"""
Forward read is a hybrid read with BC at the beginning, then linker sequence and the the insert/oligo.
This sankefile will extract the BC and FW read from the hybrid read.
"""


rule assignment_hybridFWRead_get_reads_by_length:
    """
    Get the barcode and read from the FW read using fixed length
    """
    conda:
        "../../envs/default.yaml"
    input:
        fastq=lambda wc: config["assignments"][wc.assignment]["FW"],
        check="results/assignment/{assignment}/design_check.done",
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
        zcat {input.fastq} | \
        awk '{{if (NR%4==2 || NR%4==0){{
                print substr($0,1,20) > "{output.BC_tmp}"; print substr($0,{params.insert_start}) > "{output.FW_tmp}"
            }} else {{
                print $0 > "{output.BC_tmp}"; print $0 > "{output.FW_tmp}"
            }}}}';
        cat {output.BC_tmp} | bgzip > {output.BC} & cat {output.FW_tmp} | bgzip > {output.FW};
        """


rule assignmemt_hybridFWRead_get_reads_by_cutadapt:
    """
    Get the barcode and read from the FW read using cutadapt.
    Uses the paired end mode of cutadapt to write the FW and BC read.
    """
    conda:
        "../../envs/cutadapt.yaml"
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
