######################################
### Everything before assigning BC ###
######################################

### START COUNTING ####


rule counts_onlyFW_raw_counts_by_length:
    """
    Getting the BCs from the reads using fixed length.
    """
    conda:
        "../../envs/default.yaml"
    input:
        lambda wc: getFW(wc.project, wc.condition, wc.replicate, wc.type),
    output:
        "results/experiments/{project}/counts/onlyFWByLength.{condition}_{replicate}_{type}_raw_counts.tsv.gz",
    params:
        bc_length=lambda wc: config["experiments"][wc.project]["bc_length"],
    log:
        temp(
            "results/logs/counts/onlyFW/onlyFW_raw_counts_by_length.{project}.{condition}.{replicate}.{type}.log"
        ),
    shell:
        """
        zcat {input} | \
        awk 'NR%4==2 {{print substr($1,1,{params.bc_length})}}' | \
        sort | \
        gzip -c > {output} 2> {log}
        """


rule counts_onlyFW_raw_counts_by_cutadapt:
    """
    Getting the BCs from the reads using cutadapt.
    """
    conda:
        "../../envs/cutadapt.yaml"
    input:
        lambda wc: getFW(wc.project, wc.condition, wc.replicate, wc.type),
    output:
        "results/experiments/{project}/counts/onlyFWByCutadapt.{condition}_{replicate}_{type}_raw_counts.tsv.gz",
    params:
        adapter=lambda wc: config["experiments"][wc.project]["adapter"],
    log:
        temp(
            "results/logs/counts/onlyFW/onlyFW_raw_counts_by_cutadapt.{project}.{condition}.{replicate}.{type}.log"
        ),
    shell:
        """
        zcat {input} | \
        cutadapt -a {params.adapter} - |
        awk 'NR%4==2 {{print $1}}' | \
        sort | \
        gzip -c > {output} 2> {log}
        """
