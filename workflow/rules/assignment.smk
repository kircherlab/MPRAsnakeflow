# Assignment workflow


include: "assignment/statistic.smk"
include: "assignment/assignment_common.smk"


rule assignment_get_reads_by_length:
    """
    Get the barcode and read from the FW read using fixed length
    """
    conda:
        "../envs/default.yaml"
    input:
        lambda wc: config["assignments"][wc.assignment]["FW"],
    output:
        FW_tmp=temp("results/assignment/{assignment}/fastq/FW.byLength.fastq"),
        BC_tmp=temp("results/assignment/{assignment}/fastq/BC.byLength.fastq"),
        FW="results/assignment/{assignment}/fastq/FW.byLength.fastq.gz",
        BC="results/assignment/{assignment}/fastq/BC.byLength.fastq.gz",
    log:
        temp("results/logs/assignment/get_BC_read_by_length.{assignment}.log"),
    params:
        bc_length=lambda wc: config["assignments"][wc.assignment]["bc_length"],
        insert_start=lambda wc: config["assignments"][wc.assignment]["bc_length"]
        + config["assignments"][wc.assignment]["linker_length"]
        + 1,
    shell:
        """
        zcat {input} | \
        awk '{{if (NR%4==2 || NR%4==0){{
                print substr($0,1,20) > "{output.BC_tmp}"; print substr($0,{params.insert_start}) > "{output.FW_tmp}"
            }} else {{
                print $0 > "{output.BC_tmp}"; print $0 > "{output.FW_tmp}"
            }}}}';
        cat {output.BC_tmp} | bgzip > {output.BC} & cat {output.FW_tmp} | bgzip > {output.FW};
        """


# rule assignmemt_get_read_by_cutadapt:
#     """
#     Get the barcode and read from the FW read using cutadapt
#     """
#     conda:
#         "../envs/cutadapt.yaml"
#     input:
#         lambda wc: config["assignments"][wc.assignment]["FW"],
#     output:
#         FW="results/assignment/{assignment}/fastq/FW.byCutadapt.fastq.gz",
#     log:
#         temp("results/logs/assignment/get_FW_read_by_cutadapt.{assignment}.log"),
#     params:
#         linker=lambda wc: config["assignments"][wc.assignment]["linker"],
#     shell:
#         """
#         zcat {input} | \
#         cutadapt -g {params.linker} \
#         -o {output.FW} - &> {log}
#         """


rule assignmemt_get_reads_by_cutadapt:
    """
    Get the barcode and read from the FW read using cutadapt.
    Uses the paired end mode of cutadapt to write the FW and BC read.
    """
    conda:
        "../envs/cutadapt.yaml"
    input:
        lambda wc: config["assignments"][wc.assignment]["FW"],
    output:
        BC="results/assignment/{assignment}/fastq/BC.byCutadapt.fastq.gz",
        FW="results/assignment/{assignment}/fastq/FW.byCutadapt.fastq.gz",
    log:
        temp("results/logs/assignment/get_reads_by_cutadapt.{assignment}.log"),
    params:
        linker=lambda wc: config["assignments"][wc.assignment]["linker"],
    shell:
        """
        cutadapt -a {params.linker} -G {params.linker}\
        -o {output.BC} -p {output.FW} <(zcat {input}) <(zcat {input}) &> {log}
        """


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
        BC_rev_comp= lambda wc: "--reverse-complement" if config["assignments"][wc.assignment]["BC_rev_comp"] else ""
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
        temp("results/logs/assignment/merge.{assignment}.{split}.log"),
    shell:
        """
        NGmerge \
        -1 {input.FW} \
        -2 {input.REV} \
        -m {params.min_overlap} -p {params.frac_mismatches_allowed} -e {params.min_dovetailed_overlap} \
        -z \
        -o  {output.join} \
        -i -f {output.un} \
        -l {log}
        """


assignment_bwa_dicts = ["bwt", "sa", "pac", "ann", "amb"]


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


rule assignment_mapping:
    """
    Map the reads to the reference and sort.
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


rule assignment_getBCs:
    """
    Get the barcodes.
    """
    input:
        "results/assignment/{assignment}/aligned_merged_reads.bam",
    output:
        "results/assignment/{assignment}/barcodes_incl_other.sorted.tsv.gz",
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
    log:
        temp("results/logs/assignment/getBCs.{assignment}.log"),
    shell:
        """
        samtools view -F 1792 {input} | \
        awk -v "OFS=\\t" '{{
            split($(NF),a,":");
            split(a[3],a,",");
            if (a[1] !~ /N/) {{
                if (($5 > 0) && ($4 >= {params.alignment_start_min}) && ($4 <= {params.alignment_start_max}) && (length($10) >= {params.sequence_length_min}) && (length($10) <= {params.sequence_length_max})) {{
                    print a[1],$3,$4";"$6";"$12";"$13";"$5 
                }} else {{
                    print a[1],"other","NA" 
                }}
            }}
        }}' | sort -k1,1 -k2,2 -k3,3 | gzip -c > {output} 2> {log}
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
