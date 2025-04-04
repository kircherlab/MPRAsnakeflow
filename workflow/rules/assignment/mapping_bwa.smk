rule assignment_mapping_bwa_ref:
    """
    Create mapping reference for BWA from design file.
    """
    conda:
        getCondaEnv("bwa_samtools_picard_htslib.yaml")
    input:
        ref="results/assignment/{assignment}/reference/reference.fa",
        check="results/assignment/{assignment}/design_check.done",
    output:
        bwa=expand(
            "results/assignment/{{assignment}}/reference/reference.fa.{ext}",
            ext=["fai"] + assignment_bwa_dicts,
        ),
        d="results/assignment/{assignment}/reference/reference.fa.dict",
    log:
        temp("results/logs/assignment/mapping.bwa_ref.{assignment}.log"),
    shell:
        """
        bwa index -a bwtsw {input.ref} &> {log};
        samtools faidx {input.ref} &>> {log};
        picard CreateSequenceDictionary REFERENCE={input.ref} OUTPUT={output.d} &>> {log}
        """


rule assignment_mapping_bwa:
    """
    Map the reads to the reference and sort unsing bwa mem
    """
    conda:
        getCondaEnv("bwa_samtools_picard_htslib.yaml")
    threads: 1
    input:
        reads="results/assignment/{assignment}/fastq/merge_split{split}.join.fastq.gz",
        reference="results/assignment/{assignment}/reference/reference.fa",
        bwa_index=expand(
            "results/assignment/{{assignment}}/reference/reference.fa.{ext}",
            ext=["fai", "dict"] + assignment_bwa_dicts,
        ),
    output:
        bam=temp("results/assignment/{assignment}/bwa/merge_split{split}.mapped.bam"),
    log:
        temp("results/logs/assignment/mapping.bwa.{assignment}.{split}.log"),
    shell:
        """
        bwa mem -t {threads} -L 80 -M -C {input.reference} <(
            gzip -dc {input.reads}
        )  | samtools sort -l 0 -@ {threads} > {output} 2> {log}
        """


rule assignment_mapping_bwa_getBCs:
    """
    Get the barcodes.
    """
    conda:
        getCondaEnv("bwa_samtools_picard_htslib.yaml")
    input:
        "results/assignment/{assignment}/bwa/merge_split{split}.mapped.bam",
    output:
        temp("results/assignment/{assignment}/BCs/barcodes_bwa.{split}.tsv"),
    params:
        alignment_start_min=lambda wc: config["assignments"][wc.assignment][
            "alignment_tool"
        ]["configs"]["alignment_start"]["min"],
        alignment_start_max=lambda wc: config["assignments"][wc.assignment][
            "alignment_tool"
        ]["configs"]["alignment_start"]["max"],
        sequence_length_min=lambda wc: config["assignments"][wc.assignment][
            "alignment_tool"
        ]["configs"]["sequence_length"]["min"],
        sequence_length_max=lambda wc: config["assignments"][wc.assignment][
            "alignment_tool"
        ]["configs"]["sequence_length"]["max"],
        mapping_quality_min=lambda wc: config["assignments"][wc.assignment][
            "alignment_tool"
        ]["configs"]["min_mapping_quality"],
    log:
        temp("results/logs/assignment/mapping.bwa.getBCs.{assignment}.{split}.log"),
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
    conda:
        getCondaEnv("bwa_samtools_picard_htslib.yaml")
    threads: 1
    input:
        bams=lambda wc: expand(
            "results/assignment/{{assignment}}/{mapper}/merge_split{split}.mapped.bam",
            split=range(0, getSplitNumber()),
            mapper=config["assignments"][wc.assignment]["alignment_tool"]["tool"],
        ),
    output:
        "results/assignment/{assignment}/aligned_merged_reads.bam",
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
    conda:
        getCondaEnv("bwa_samtools_picard_htslib.yaml")
    input:
        "results/assignment/{assignment}/aligned_merged_reads.bam",
    output:
        "results/assignment/{assignment}/aligned_merged_reads.bam.bai",
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
    conda:
        getCondaEnv("bwa_samtools_picard_htslib.yaml")
    input:
        bam="results/assignment/{assignment}/aligned_merged_reads.bam",
        idx="results/assignment/{assignment}/aligned_merged_reads.bam.bai",
    output:
        "results/assignment/{assignment}/statistic/assignment/bam_stats.txt",
    log:
        temp("results/logs/assignment/flagstat.{assignment}.log"),
    shell:
        """
        samtools flagstat {input.bam} > {output} 2> {log}
        """
