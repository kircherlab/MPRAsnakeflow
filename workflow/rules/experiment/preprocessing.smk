rule experiment_preprocessing_split_reads:
    """
Split the fastq files into n files for parallelisation.
n is given by split_read in the configuration file.
"""
    input:
        lambda wc: getExperimentReads(
            wc.read_type, wc.project, wc.condition, wc.replicate, wc.type, check_splitting=False, check_trimming=False
        ),
    output:
        temp(
            expand(
                "results/experiments/{{project}}/fastq/{{read_type}}.split.{{condition}}.{{replicate}}.{{type}}.{split}.fastq.gz",
                split=range(
                    0,
                    getMaxExperimentSplitNumber(),
                ),
            ),
        ),
    log:
        "results/logs/experiment/preprocessing/split_reads.{project}.{condition}.{replicate}.{type}.{read_type}.log",
    conda:
        getCondaEnv("fastqsplitter.yaml")
    shell:
        """
        OUTPUT_FILES=""
        for file in {output}; do
            OUTPUT_FILES="${{OUTPUT_FILES}} -o ${{file}}"
        done
        fastqsplitter -i <(zcat {input}) -t 1 ${{OUTPUT_FILES}} &> {log}
        """


rule experiment_preprocessing_trim_reads:
    """
Getting the BCs from the reads using cutadapt.
"""
    input:
        lambda wc: getExperimentReads(
            wc.read_type, wc.project, wc.condition, wc.replicate, wc.type, check_splitting=True, check_trimming=False
        ),
    output:
        trimmed_reads="results/experiments/{project}/fastq/{read_type}.trimmed.{condition}.{replicate}.{type}.{split}.fastq.gz",
    log:
        "results/logs/experiment/preprocessing/trim_reads.{project}.{condition}.{replicate}.{type}.{split}.{read_type}.log",
    wildcard_constraints:
        read_type=r"(FWD)|(REV)|(UMI)",
    conda:
        getCondaEnv("cutadapt.yaml")
    threads: 1
    params:
        adapters=lambda wc: getExperimentCutadaptAdapters(wc.project, wc.read_type),
    shell:
        """
        cutadapt --cores {threads} {params.adapters} \
        -o {output.trimmed_reads} <(zcat {input.reads}) &> {log}
        """
