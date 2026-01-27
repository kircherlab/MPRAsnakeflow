"""
Forward read is a hybrid read with BC at the beginning, then linker sequence and the the insert/oligo.
This sankefile will extract the BC and FW read from the hybrid read.
"""


rule assignment_hybridFWDRead_get_reads_by_length:
    """
    Get the barcode and read from the FW read using fixed length
    """
    conda:
        getCondaEnv("default.yaml")
    input:
        fastq=lambda wc: (
            "results/assignment/{assignment}/fastq/FWD.trimmed.fastq.gz"
            if useAssignmentAdapterTrimming(wc.assignment, "FWD")
            else config["assignments"][wc.assignment]["FWD"]
        ),
        check="results/assignment/{assignment}/design_check.done",
    output:
        FWD_tmp=temp("results/assignment/{assignment}/fastq/FWD.byLength.fastq"),
        BC_tmp=temp("results/assignment/{assignment}/fastq/BC.byLength.fastq"),
        FWD="results/assignment/{assignment}/fastq/FWD.byLength.fastq.gz",
        BC="results/assignment/{assignment}/fastq/BC.byLength.fastq.gz",
    log:
        temp(
            "results/logs/assignment/hybridFWDRead_get_reads_by_length.{assignment}.log"
        ),
    params:
        bc_length=lambda wc: config["assignments"][wc.assignment]["bc_length"],
        insert_start=lambda wc: config["assignments"][wc.assignment]["bc_length"]
        + config["assignments"][wc.assignment]["linker_length"]
        + 1,
    shell:
        """
        zcat {input.fastq} | \
        awk '{{if (NR%4==2 || NR%4==0){{
                print substr($0,1,{params.bc_length}) > "{output.BC_tmp}"; print substr($0,{params.insert_start}) > "{output.FWD_tmp}"
            }} else {{
                print $0 > "{output.BC_tmp}"; print $0 > "{output.FWD_tmp}"
            }}}}';
        cat {output.BC_tmp} | bgzip > {output.BC} & cat {output.FWD_tmp} | bgzip > {output.FWD};
        """


rule assignment_hybridFWDRead_get_reads_by_cutadapt:
    """
    Get the barcode and read from the FWD read using cutadapt.
    Uses the paired end mode of cutadapt to write the FWD and BC read.
    """
    conda:
        getCondaEnv("cutadapt.yaml")
    threads: 1
    input:
        lambda wc: (
            "results/assignment/{assignment}/fastq/FWD.trimmed.fastq.gz"
            if useAssignmentAdapterTrimming(wc.assignment, "FWD")
            else config["assignments"][wc.assignment]["FWD"]
        ),
    output:
        BC="results/assignment/{assignment}/fastq/BC.byCutadapt.fastq.gz",
        FWD="results/assignment/{assignment}/fastq/FWD.byCutadapt.fastq.gz",
    log:
        temp(
            "results/logs/assignment/hybridFWDRead_get_reads_by_cutadapt.{assignment}.log"
        ),
    params:
        linker=lambda wc: config["assignments"][wc.assignment]["linker"],
    shell:
        """
        cutadapt --cores {threads} -a {params.linker} -G {params.linker}\
        -o {output.BC} -p {output.FWD} <(zcat {input}) <(zcat {input}) &> {log}
        """
