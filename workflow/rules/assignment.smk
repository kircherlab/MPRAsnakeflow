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
        un="results/assignment/{assignment}/fastq/merge_split{split}.un.fastq.gz",
        join="results/assignment/{assignment}/fastq/merge_split{split}.join.fastq.gz",
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
        lambda wc: config["assignments"][wc.assignment]["design_file"],
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
    log:
        temp("results/logs/assignment/mapping.{assignment}.{split}.log"),
    shell:
        """
        bwa mem -t {threads} -L 80 -M -C {input.reference} <(
            gzip -dc {input.reads}
        )  | samtools sort -l 0 -@ {threads} > {output} 2> {log}
        """

rule assignment_getBCs_additional_filter:
    """
    Get the barcodes with a python script to rescue alignments with 0 mapping quality according to bwa.
    """
    input:
        bam="results/assignment/{assignment}/bam/merge_split{split}.mapped.bam",
        script=getScript("assignment/filter_bc_from_bam.py"),
    output:
        "results/assignment/{assignment}/BCs/barcodes_incl_other.{split}.tsv",
    params:
        identity_threshold=0.98,
        mismatches_threshold=3,
        use_expected_alignment_length=True,
        expected_alignment_length=265,
        verbose=True,
    conda:
        "../envs/getBC.yaml"
    log:
        "results/logs/assignment/getBCs.{assignment}.{split}.log",
    shell:
        """
        python {input.script} \
        --identity_threshold {params.identity_threshold} --mismatches_threshold {params.mismatches_threshold} \
        --use_expected_alignment_length {params.use_expected_alignment_length} --expected_alignment_length {params.expected_alignment_length} \
        --bamfile {input.bam} --verbose {params.verbose} --output {output} 2> {log} && sort -k1,1 -k2,2 -k3,3 -o {output} {output}
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
        expand(
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
        sort --batch-size={params.batch_size} --parallel={threads} -k1,1 -k2,2 -k3,3 -m {input} | gzip -c > {output} 2> {log}
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
        unknown_other="-o",
        ambiguous="-a",
        bc_length=lambda wc: config["assignments"][wc.assignment]["bc_length"],
    shell:
        """
        zcat  {input.assignment} | \
        awk -v "OFS=\\t" -F"\\t" '{{if (length($1)=={params.bc_length}){{print $0 }}}}' | \
        python {input.script} \
        -m {params.min_support} -f {params.fraction} {params.unknown_other} {params.ambiguous} | \
        gzip -c > {output} 2> {log}
        """
