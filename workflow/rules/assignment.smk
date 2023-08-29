"""
Assignment workflow

This workflow will asssign barcodes to the designed reference/insert/oligos.
The output is a tabular file that matched barcodes with oligos.
"""


include: "assignment/common.smk"
include: "assignment/hybridFWRead.smk"
include: "assignment/statistic.smk"


rule assignment_fastq_split:
    """
    Split the fastq files into n files for parallelisation. 
    n is given by split_read in the configuration file.
    """
    input:
        lambda wc: getAssignmentRead(wc.assignment, wc.read),
    output:
        temp(
            expand(
                "results/assignment/{{assignment}}/fastq/splits/{{read}}.split{split}.fastq.gz",
                split=range(0, getSplitNumber()),
            ),
        ),
    conda:
        "../envs/fastqsplitter.yaml"
    log:
        temp("results/logs/assignment/fastq_split.{assignment}.{read}.log"),
    params:
        files=lambda wc: " ".join(
            [
                "-o %s" % i
                for i in expand(
                    "results/assignment/{assignment}/fastq/splits/{read}.split{split}.fastq.gz",
                    assignment=wc.assignment,
                    read=wc.read,
                    split=range(0, getSplitNumber()),
                )
            ]
        ),
    shell:
        """
        fastqsplitter -i <(zcat {input}) -t 1 {params.files} &> {log}
        """


rule assignment_attach_idx:
    """
    Extract the index sequence and add it to the header.
    """
    conda:
        "../envs/NGmerge.yaml"
    input:
        read="results/assignment/{assignment}/fastq/splits/{read}.split{split}.fastq.gz",
        BC="results/assignment/{assignment}/fastq/splits/BC.split{split}.fastq.gz",
        script=getScript("attachBCToFastQ.py"),
    output:
        read=temp(
            "results/assignment/{assignment}/fastq/splits/{read}.split{split}.BCattached.fastq.gz"
        ),
    params:
        BC_rev_comp=lambda wc: "--reverse-complement"
        if config["assignments"][wc.assignment]["BC_rev_comp"]
        else "",
    log:
        temp("results/logs/assignment/attach_idx.{assignment}.{split}.{read}.log"),
    shell:
        """
        python {input.script} -r {input.read} -b {input.BC} {params.BC_rev_comp} | bgzip -c > {output.read} 2> {log}
        """


rule assignment_merge:
    """
    Merge the FW,REV and BC fastq files into one. 
    Extract the index sequence and add it to the header.
    """
    conda:
        "../envs/NGmerge.yaml"
    input:
        FW="results/assignment/{assignment}/fastq/splits/FW.split{split}.BCattached.fastq.gz",
        REV="results/assignment/{assignment}/fastq/splits/REV.split{split}.BCattached.fastq.gz",
    output:
        un=temp("results/assignment/{assignment}/fastq/merge_split{split}.un.fastq.gz"),
        join=temp(
            "results/assignment/{assignment}/fastq/merge_split{split}.join.fastq.gz"
        ),
    params:
        min_overlap=lambda wc: config["assignments"][wc.assignment]["NGmerge"][
            "min_overlap"
        ],
        frac_mismatches_allowed=lambda wc: config["assignments"][wc.assignment][
            "NGmerge"
        ]["frac_mismatches_allowed"],
        min_dovetailed_overlap=lambda wc: config["assignments"][wc.assignment][
            "NGmerge"
        ]["min_dovetailed_overlap"],
    log:
        temp("results/logs/assignment/merge.{assignment}.{split}.log.gz"),
    shell:
        """
        NGmerge \
        -1 {input.FW} \
        -2 {input.REV} \
        -m {params.min_overlap} -p {params.frac_mismatches_allowed} \
        -d \
        -e {params.min_dovetailed_overlap} \
        -z \
        -o  {output.join} \
        -i -f {output.un} \
        -l >(gzip -c - > {log})
        """


rule assignment_bwa_ref:
    """
    Create mapping reference for BWA from design file.
    """
    input:
        lambda wc: config["assignments"][wc.assignment]["reference"],
    output:
        ref="results/assignment/{assignment}/reference/reference.fa",
        bwa=expand(
            "results/assignment/{{assignment}}/reference/reference.fa.{ext}",
            ext=["fai"] + assignment_bwa_dicts,
        ),
        d="results/assignment/{assignment}/reference/reference.fa.dict",
    conda:
        "../envs/bwa_samtools_picard_htslib.yaml"
    log:
        temp("results/logs/assignment/bwa_ref.{assignment}.log"),
    shell:
        """
        cat {input} | awk '{{gsub(/[\\]\\[]/,"_")}}$0' > {output.ref};
        bwa index -a bwtsw {output.ref} &> {log};
        samtools faidx {output.ref} &>> {log};
        picard CreateSequenceDictionary REFERENCE={output.ref} OUTPUT={output.d} &>> {log}
        """


rule assignment_mapping_bwa:
    """
    Map the reads to the reference and sort unsing bwa mem
    """
    input:
        reads="results/assignment/{assignment}/fastq/merge_split{split}.join.fastq.gz",
        reference="results/assignment/{assignment}/reference/reference.fa",
        bwa_index=expand(
            "results/assignment/{{assignment}}/reference/reference.fa.{ext}",
            ext=["fai", "dict"] + assignment_bwa_dicts,
        ),
    output:
        bam=temp("results/assignment/{assignment}/bam/merge_split{split}.mapped.bam"),
    conda:
        "../envs/bwa_samtools_picard_htslib.yaml"
    threads: config["global"]["threads"]
    log:
        temp("results/logs/assignment/mapping.{assignment}.{split}.log"),
    shell:
        """
        bwa mem -t {threads} -L 80 -M -C {input.reference} <(
            gzip -dc {input.reads}
        )  | samtools sort -l 0 -@ {threads} > {output} 2> {log}
        """


rule assignment_getBCs:
    """
    Get the barcodes.
    """
    input:
        "results/assignment/{assignment}/bam/merge_split{split}.mapped.bam",
    output:
        temp("results/assignment/{assignment}/BCs/barcodes_incl_other.{split}.tsv"),
    conda:
        "../envs/bwa_samtools_picard_htslib.yaml"
    params:
        alignment_start_min=lambda wc: config["assignments"][wc.assignment][
            "alignment_start"
        ]["min"],
        alignment_start_max=lambda wc: config["assignments"][wc.assignment][
            "alignment_start"
        ]["max"],
        sequence_length_min=lambda wc: config["assignments"][wc.assignment][
            "sequence_length"
        ]["min"],
        sequence_length_max=lambda wc: config["assignments"][wc.assignment][
            "sequence_length"
        ]["max"],
        mapping_quality_min=lambda wc: config["assignments"][wc.assignment][
            "min_mapping_quality"
        ],
    log:
        temp("results/logs/assignment/getBCs.{assignment}.{split}.log"),
    shell:
        """
        export LC_ALL=C # speed up sorting
        samtools view -F 1792 {input} | \
        awk -v "OFS=\\t" '{{
            split($(NF),a,":");
            split(a[3],a,",");
            if (a[1] !~ /N/) {{
                if (($5 >= {params.mapping_quality_min}) && ($4 >= {params.alignment_start_min}) && ($4 <= {params.alignment_start_max}) && (length($10) >= {params.sequence_length_min}) && (length($10) <= {params.sequence_length_max})) {{
                    print a[1],$3,$4";"$6";"$12";"$13";"$5 
                }} else {{
                    print a[1],"other","NA" 
                }}
            }}
        }}' | sort -k1,1 -k2,2 -k3,3 -S 7G > {output} 2> {log}
        """


rule assignment_collect:
    """
    Collect mapped reads.
    """
    input:
        bams=expand(
            "results/assignment/{{assignment}}/bam/merge_split{split}.mapped.bam",
            split=range(0, getSplitNumber()),
        ),
    output:
        "results/assignment/{assignment}/aligned_merged_reads.bam",
    conda:
        "../envs/bwa_samtools_picard_htslib.yaml"
    threads: config["global"]["threads"]
    log:
        temp("results/logs/assignment/collect.{assignment}.log"),
    shell:
        """
        samtools merge -@ {threads} {output} {input.bams} 2> {log}
        """


rule assignment_idx_bam:
    """
    Index the BAM file
    """
    input:
        "results/assignment/{assignment}/aligned_merged_reads.bam",
    output:
        "results/assignment/{assignment}/aligned_merged_reads.bam.bai",
    conda:
        "../envs/bwa_samtools_picard_htslib.yaml"
    log:
        "results/logs/assignment/{assignment}/assignment_idx_bam.log",
    shell:
        """
        samtools index {input} 2> {log}
        """


rule assignment_flagstat:
    """
    Run samtools flagstat
    """
    input:
        bam="results/assignment/{assignment}/aligned_merged_reads.bam",
        idx="results/assignment/{assignment}/aligned_merged_reads.bam.bai",
    output:
        "results/assignment/{assignment}/statistic/assignment/bam_stats.txt",
    conda:
        "../envs/bwa_samtools_picard_htslib.yaml"
    log:
        temp("results/logs/assignment/flagstat.{assignment}.log"),
    shell:
        """
        samtools flagstat {input.bam} > {output} 2> {log}
        """


rule assignment_collectBCs:
    """
    Get the barcodes.
    """
    input:
        lambda wc: expand(
            "results/assignment/{{assignment}}/BCs/barcodes_exact.{split}.tsv",
                split=range(0, getSplitNumber()),
            )
            if config["assignments"][wc.assignment]["alignment_tool"]["tool"] == "exact"
        else expand(
            "results/assignment/{{assignment}}/BCs/barcodes_incl_other.{split}.tsv",
            split=range(0, getSplitNumber()),
        ),
    output:
        "results/assignment/{assignment}/barcodes_incl_other.sorted.tsv.gz",
    params:
        batch_size=getSplitNumber(),
    threads: 20
    conda:
        "../envs/default.yaml"
    log:
        temp("results/logs/assignment/collectBCs.{assignment}.log"),
    shell:
        """
        export LC_ALL=C # speed up sort
        sort -S 7G --batch-size={params.batch_size} --parallel={threads} -k1,1 -k2,2 -k3,3 -m {input} | \
        gzip -c > {output} 2> {log}
        """


rule assignment_filter:
    """
    Filter the barcodes file based on the config given in the config-file.
    """
    input:
        assignment="results/assignment/{assignment}/barcodes_incl_other.sorted.tsv.gz",
        script=getScript("assignment/filterAssignmentTsv.py"),
    output:
        "results/assignment/{assignment}/assignment_barcodes.{assignment_config}.sorted.tsv.gz",
    conda:
        "../envs/python3.yaml"
    log:
        temp("results/logs/assignment/filter.{assignment}.{assignment_config}.log"),
    params:
        min_support=lambda wc: config["assignments"][wc.assignment]["configs"][
            wc.assignment_config
        ]["min_support"],
        fraction=lambda wc: config["assignments"][wc.assignment]["configs"][
            wc.assignment_config
        ]["fraction"],
        unknown_other=lambda wc: "-o"
        if config["assignments"][wc.assignment]["configs"][wc.assignment_config][
            "unknown_other"
        ]
        else "",
        ambiguous=lambda wc: "-a"
        if config["assignments"][wc.assignment]["configs"][wc.assignment_config][
            "ambiguous"
        ]
        else "",
        bc_length=lambda wc: config["assignments"][wc.assignment]["bc_length"],
    shell:
        """
        zcat  {input.assignment} | \
        awk -v "OFS=\\t" -F"\\t" '{{if (length($1)=={params.bc_length}){{print $0 }}}}' | \
        python {input.script} \
        -m {params.min_support} -f {params.fraction} {params.unknown_other} {params.ambiguous} | \
        gzip -c > {output} 2> {log}
        """

rule assignment_mapping_exact_reference:
    """
    Create reference to map the exact design
    """
    input:
        lambda wc: config["assignments"][wc.assignment]["reference"],
    output:
        "results/assignment/{assignment}/reference/reference_exact.fa",
    shell:
        """
        paste <(
            cat {input} | awk '{{if ($1 ~ /^>/) {{ gsub(/[\\]\\[]/,"_"); print substr($1,2)}}}}';
            cat {input} | awk '{{if ($1 ~ /^>/) {{ gsub(/[\\]\\[]/,"_"); print substr($1,2)}}}}';
        ) <(
            cat {input} | awk '{{if ($1 ~ /^[^>]/) {{ seq=seq$1}}; if ($1 ~ /^>/ && NR!=1) {{print seq; seq=""}}}} END {{print seq}}';
            cat {input} | awk '{{if ($1 ~ /^[^>]/) {{ seq=seq$1}}; if ($1 ~ /^>/ && NR!=1) {{print seq; seq=""}}}} END {{print seq}}' | tr ACGTacgt TGCAtgca | rev;
        ) > {output}
        """

# TODO: Set correct length using config
rule assignment_mapping_exact:
    """
    Map the reads to the reference and sort using exact match.
    """
    input:
        reads="results/assignment/{assignment}/fastq/merge_split{split}.join.fastq.gz",
        reference="results/assignment/{assignment}/reference/reference_exact.fa",
    output:
        temp("results/assignment/{assignment}/BCs/barcodes_exact.{split}.tsv"),
    conda:
        "../envs/default.yaml"
    log:
        temp("results/logs/assignment/mapping_exact.{assignment}.{split}.log"),
    shell:
        """
        # Look up exact matches in design file
        export LC_ALL=C # speed up sort

        awk -v "OFS=\\t" 'NR==FNR {{a[$2] = $1; next}} {{if ($3 in a) print $2,a[$3],"200M"; else print $2,"other","NA"}}' \
        <(
            cat {input.reference} | awk -v "OFS=\\t" '{{print $1,substr($2, 16,200)}}'
        ) \
        <(
            zcat {input.reads} | awk 'NR%4==2 || NR%4==1' | paste - -
        ) | \
        sed 's/,YI:Z[^\\t]*//g' | sed 's/XI:Z://g' | \
        sort -k1,1 -k2,2 -k3,3 -S 7G > {output} 2> {log}
        """
