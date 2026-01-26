rule assignment_preprocessing_adapter_remove:
    """
    Remove adapter sequence from the reads (3' or 5').
    Uses cutadapt to trim adapters based on the primer direction.
    """
    conda:
        getCondaEnv("cutadapt.yaml")
    threads: 1
    input:
        reads=lambda wc: config["assignments"][wc.assignment][wc.read],
    output:
        trimmed_reads=temp(
            "results/assignment/{assignment}/fastq/{read}.trimmed.fastq.gz"
        ),
    wildcard_constraints:
        read=r"(FWD)|(REV)|(BC)",
    params:
        adapters=lambda wc: getAssignmentCutadaptAdapters(wc.assignment, wc.read),
    log:
        "results/logs/assignment/preprocessing/adapter_remove.{assignment}.{read}.log",
    shell:
        """
        cutadapt --cores {threads} {params.adapters} \
        -o {output.trimmed_reads} <(zcat {input.reads}) &> {log}
        """
