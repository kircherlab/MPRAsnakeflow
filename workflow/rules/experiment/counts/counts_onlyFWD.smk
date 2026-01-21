######################################
### Everything before assigning BC ###
######################################

### START COUNTING ####


rule experiment_counts_onlyFWD_raw_counts:
    """
    Getting the BCs from the reads using fixed length.
    """
    conda:
        getCondaEnv("default.yaml")
    input:
        lambda wc: getFWD(
            wc.project, wc.condition, wc.replicate, wc.type, check_trimming=True
        ),
    output:
        "results/experiments/{project}/counts/onlyFW.{condition}.{replicate}.{type}.raw_counts.tsv.gz",
    params:
        bc_length=lambda wc: config["experiments"][wc.project]["bc_length"],
        bc_extraction=lambda wc: config["experiments"][wc.project].get(
            "bc_extraction", "start"
        ),
    log:
        temp(
            "results/logs/experiment/counts/onlyFWD/onlyFWD_raw_counts.{project}.{condition}.{replicate}.{type}.log"
        ),
    shell:
        """
        zcat {input} | \
        awk 'NR%4==2 {{if ("{params.bc_extraction}" == "start") print substr($1,1,{params.bc_length}); else print substr($1,length($1)-{params.bc_length}+1)}}' |\
        sort | \
        gzip -c > {output} 2> {log}
        """
