rule experiment_preprocessing_trim_reads:
    """
    Getting the BCs from the reads using cutadapt.
    """
    conda:
        getCondaEnv("cutadapt.yaml")
    threads: 1
    input:
        lambda wc: getReads(
            wc.read_type, wc.project, wc.condition, wc.replicate, wc.type
        ),
    output:
        trimmed_reads="results/experiments/{project}/fastq/{read_type}.trimmed.{condition}.{replicate}.{type}.fastq.gz",
    wildcard_constraints:
        read_type=r"(FWD)|(REV)|(UMI)",
    params:
        adapters=lambda wc: getExperimentCutadaptAdapters(wc.project, wc.read_type),
    log:
        "results/logs/experiment/preprocessing/trim_reads.{project}.{condition}.{replicate}.{type}.{read_type}.log",
    shell:
        """
        cutadapt --cores {threads} {params.adapters} \
        -o {output.trimmed_reads} <(zcat {input.reads}) &> {log}
        """
