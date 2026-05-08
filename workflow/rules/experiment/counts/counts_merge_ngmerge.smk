######################################
### Shared merge rule templates    ###
######################################


rule experiment_counts_merge_NGmerge_template:
    """
Template rule to merge paired reads using NGmerge.
"""
    log:
        temp("results/logs/experiment/counts/merge_NGmerge.template.{project}.{condition}.{replicate}.{type}.{split}.log"),
    conda:
        getCondaEnv("NGmerge.yaml")
    params:
        min_overlap=lambda wc: config["experiments"][wc.project].get("NGmerge", {}).get("min_overlap", 11),
        frac_mismatches_allowed=lambda wc: config["experiments"][wc.project]
        .get("NGmerge", {})
        .get("frac_mismatches_allowed", 0.1),
    shell:
        """
        NGmerge \
        -1 {input.FWD} \
        -2 {input.REV} \
        -m {params.min_overlap} -p {params.frac_mismatches_allowed} \
        -z \
        -o  {output.join} \
        -i -f {output.un} &> {log}
        """
