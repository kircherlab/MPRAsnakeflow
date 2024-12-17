######################################
### Everything before assigning BC ###
######################################


### Create_BAM_umi without demultiplexing ###


rule counts_onlyFWUMI_raw_counts:
    """
    Getting the BCs and UMIs from the reads using fixed length.
    """
    conda:
        "../../envs/default.yaml"
    input:
        fw_fastq=lambda wc: getFW(wc.project, wc.condition, wc.replicate, wc.type),
        umi_fastq=lambda wc: getUMI(wc.project, wc.condition, wc.replicate, wc.type),
    output:
        "results/experiments/{project}/counts/onlyFWUMI.{condition}_{replicate}_{type}_raw_counts.tsv.gz",
    params:
        umi_length=lambda wc: config["experiments"][wc.project]["umi_length"],
    log:
        temp(
            "results/logs/counts/onlyFW/onlyFWUMI_raw_counts_by_length.{project}.{condition}.{replicate}.{type}.log"
        ),
    shell:
        """
        paste <(zcat {input.fw_fastq} | awk 'NR%4==2 {{print $1}}') \
              <(zcat {input.umi_fastq} | awk 'NR%4==2 {{print $1}}') | \
        awk -v 'OFS=\\t' 'length($2) == {params.umi_length} {{print $0}}' | \
        sort | uniq -c | \
        awk -v 'OFS=\\t' '{{ print $2,$3,$1 }}' | \
        gzip -c > {output} 2> {log}
        """
