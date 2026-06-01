rule assignment_mapping_pbmm2_index:
    """
Create pbmm2 index from design reference.
"""
    input:
        ref="results/assignment/{assignment}/reference/reference.fa",
        check="results/assignment/{assignment}/design_check.done",
    output:
        "results/assignment/{assignment}/reference/reference.fa.mmi",
    log:
        temp("results/logs/assignment/mapping_pbmm2_index.{assignment}.log"),
    conda:
        getCondaEnv("pbmm2_pysam.yaml")
    shell:
        """
        pbmm2 index {input.ref} {output} &>{log}
        """


rule assignment_mapping_pbmm2_align:
    """
Align long reads (BAM or FASTA) to reference using pbmm2.
"""
    input:
        reads=lambda wc: config["assignments"][wc.assignment]["long_read_input"],
        index="results/assignment/{assignment}/reference/reference.fa.mmi",
    output:
        "results/assignment/{assignment}/pbmm2/aligned.bam",
    log:
        temp("results/logs/assignment/mapping_pbmm2_align.{assignment}.log"),
    conda:
        getCondaEnv("pbmm2_pysam.yaml")
    threads: 8
    params:
        preset=lambda wc: config["assignments"][wc.assignment]["alignment_tool"]["configs"]["preset"],
        min_concordance=lambda wc: config["assignments"][wc.assignment]["alignment_tool"]["configs"]["min_concordance"] * 100,
    shell:
        """
        pbmm2 align {input.index} {input.reads} {output} \
            --preset {params.preset} \
            --sort \
            --best-n 1 \
            --min-concordance-perc {params.min_concordance} \
            -j {threads} &>{log}
        """


rule assignment_mapping_pbmm2_getBCs:
    """
Extract barcodes from aligned long reads. Produces the standard
barcode TSV for downstream collection and filtering.
"""
    input:
        bam="results/assignment/{assignment}/pbmm2/aligned.bam",
        script=getScript("assignment/longread_extract.py"),
    output:
        temp("results/assignment/{assignment}/BCs/barcodes.pbmm2.tsv"),
    log:
        temp("results/logs/assignment/mapping_pbmm2_getBCs.{assignment}.log"),
    conda:
        getCondaEnv("pbmm2_pysam.yaml")
    params:
        pattern=lambda wc: config["assignments"][wc.assignment]["linker"],
        bc_length=lambda wc: config["assignments"][wc.assignment]["bc_length"],
    shell:
        """
        export LC_ALL=C
        python {input.script} \
            --bam {input.bam} \
            --pattern {params.pattern} \
            --bc-length {params.bc_length} \
            --output /dev/stdout 2>{log} \
            | sort -k1,1 -k2,2 -k3,3 -S 7G >{output}
        """
