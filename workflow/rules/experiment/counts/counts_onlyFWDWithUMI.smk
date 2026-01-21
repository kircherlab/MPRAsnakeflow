######################################
### Everything before assigning BC ###
######################################


### Create_BAM_umi without demultiplexing ###


rule experiment_counts_onlyFWDUMI_raw_counts:
    """
    Getting the BCs and UMIs from the reads using fixed length.
    """
    conda:
        getCondaEnv("default.yaml")
    input:
        fw_fastq=lambda wc: getFWD(
            wc.project, wc.condition, wc.replicate, wc.type, check_trimming=True
        ),
        umi_fastq=lambda wc: getUMI(
            wc.project, wc.condition, wc.replicate, wc.type, check_trimming=True
        ),
    output:
        "results/experiments/{project}/counts/onlyFWDUMI.{condition}.{replicate}.{type}.raw_counts.tsv.gz",
    params:
        bc_length=lambda wc: config["experiments"][wc.project]["bc_length"],
        umi_length=lambda wc: config["experiments"][wc.project]["umi_length"],
        bc_extraction=lambda wc: config["experiments"][wc.project].get(
            "bc_extraction", "start"
        ),
        umi_extraction=lambda wc: config["experiments"][wc.project].get(
            "umi_extraction", "start"
        ),
    log:
        temp(
            "results/logs/experiment/counts/onlyFWD/onlyFWDUMI_raw_counts_by_length.{project}.{condition}.{replicate}.{type}.log"
        ),
    shell:
        """
        paste <(zcat {input.fw_fastq} | awk 'NR%4==2 {{if ("{params.bc_extraction}" == "start") print substr($1,1,{params.bc_length}); else print substr($1,length($1)-{params.bc_length}+1)}}') \
              <(zcat {input.umi_fastq} | awk 'NR%4==2 {{if ("{params.umi_extraction}" == "start") print substr($1,1,{params.umi_length}); else print substr($1,length($1)-{params.umi_length}+1)}}') | \
        awk -v 'OFS=\\t' 'length($2) == {params.umi_length} {{print $0}}' | \
        sort | uniq -c | \
        awk -v 'OFS=\\t' '{{ print $2,$3,$1 }}' | \
        gzip -c > {output} 2> {log}
        """
