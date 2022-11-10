include: "assignment/statistic.smk"


rule assignment_fastq_split:
    """
    Split the fastq files into n files for parallelisation. 
    n is given by split_read in the configuration file.
    """
    input:
        lambda wc: config["assignments"][wc.assignment][wc.R],
    output:
        temp(
            expand(
                "results/assignment/{{assignment}}/fastq/splits/{{R}}.split{split}.fastq.gz",
                split=range(0, getSplitNumber()),
            ),
        ),
    conda:
        "../envs/fastqsplitter.yaml"
    log:
        temp("results/logs/assignment/fastq_split.{assignment}.{R}.log"),
    params:
        files=lambda wc: " ".join(
            [
                "-o %s" % i
                for i in expand(
                    "results/assignment/{assignment}/fastq/splits/{R}.split{split}.fastq.gz",
                    assignment=wc.assignment,
                    R=wc.R,
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
        read="results/assignment/{assignment}/fastq/splits/R{R}.split{split}.fastq.gz",
        BC="results/assignment/{assignment}/fastq/splits/R2.split{split}.fastq.gz",
        script=getScript("attachBCToFastQ.py"),
    output:
        read=temp(
            "results/assignment/{assignment}/fastq/splits/R{R}.split{split}.BCattached.fastq.gz"
        ),
    log:
        temp("results/logs/assignment/attach_idx.{assignment}.{split}.{R}.log"),
    shell:
        """
        python {input.script} -r {input.read} -b {input.BC} | bgzip -c > {output.read}
        """


rule assignment_merge:
    """
    Merge the FW,REV and BC fastq files into one. 
    Extract the index sequence and add it to the header.
    """
    conda:
        "../envs/NGmerge.yaml"
    input:
        R1="results/assignment/{assignment}/fastq/splits/R1.split{split}.BCattached.fastq.gz",
        R3="results/assignment/{assignment}/fastq/splits/R3.split{split}.BCattached.fastq.gz",
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
        -1 {input.R1} \
        -2 {input.R3} \
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
    shell:
        """
        zcat  {input.assignment} | \
        python {input.script} \
        -m {params.min_support} -f {params.fraction} {params.unknown_other} {params.ambiguous} | \
        gzip -c > {output} 2> {log}
        """
