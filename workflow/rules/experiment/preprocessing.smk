rule experiment_preprocess_trim_reads:
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
        trimmed_reads=temp(
            "results/experiments/{project}/fastq/{read_type}.trimmed.{condition}.{replicate}.{type}.fastq.gz"
        ),
    params:
        adapter=lambda wc: config["experiments"][wc.project]["adapters"],
    wildcard_constraints:
        read_type=r"(FWD)|(REV)|(UMI)",
    params:
        adapters=lambda wc: getExperimentCutadaptAdapters(wc.project, wc.read_type),
    log:
        temp(
            "results/logs/experiment/preprocess/trim_reads.{read_type}.{project}.{condition}.{replicate}.{type}.log"
        ),
    shell:
        """
        cutadapt --cores {threads} {params.adapters} \
        -o {output.trimmed_reads} <(zcat {input.reads}) &> {log}
        """
