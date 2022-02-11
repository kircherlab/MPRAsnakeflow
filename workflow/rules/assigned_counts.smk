##################################
### Assign BCs and afterwards ###
##################################


rule assignBarcodes:
    """
    Assign RNA and DNA barcodes seperately to make the statistic for assigned
    """
    conda:
        "../envs/mpraflow_py36.yaml"
    input:
        counts="results/{project}/counts/{condition}_{replicate}_{type}_final_counts_{sampling}.tsv.gz",
        association=lambda wc: config[wc.project]["assignments"][wc.assignment],
    output:
        counts="results/{project}/assigned_counts/{assignment}/{condition}_{replicate}_{type}_final_counts_{sampling}.tsv.gz",
        stats="results/{project}/stats/assigned_counts/{assignment}/{condition}_{replicate}_{type}_{sampling}.statistic.tsv.gz",
    params:
        name="{condition}_{replicate}_{type}",
    shell:
        """
        python {SCRIPTS_DIR}/count/merge_BC_and_assignment.py --counts {input.counts} \
        --assignment {input.association} \
        --output {output.counts} \
        --statistic {output.stats} \
        --name {params.name}
        """


rule createAssignmentPickleFile:
    conda:
        "envs/mpraflow_py36.yaml"
    input:
        lambda wc: config[wc.project]["assignments"][wc.assignment],
    output:
        "results/{project}/assigned_counts/{assignment}/assignment.pickle",
    shell:
        """
        python {SCRIPTS_DIR}/count/create_pickle.py -i {input} -o {output}
        """


rule dna_rna_merge:
    conda:
        "../envs/mpraflow_py36.yaml"
    input:
        counts=lambda wc: expand(
            "results/{{project}}/counts/merged/{mergeType}/{{condition}}_{{replicate}}_merged_counts_{{sampling}}.tsv.gz",
        mergeType="withoutZeros"
            if config[wc.project]["configs"][wc.config]["minRNACounts"] > 0
            and config[wc.project]["configs"][wc.config]["minDNACounts"] > 0
            else "withZeros",
        ),
        association=lambda wc: config[wc.project]["assignments"][wc.assignment],
    output:
        counts="results/{project}/assigned_counts/{assignment}/{config}/{condition}_{replicate}_merged_assigned_counts_{sampling}.tsv.gz",
        stats="results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_{replicate}_merged_assigned_counts_{sampling}.statistic.tsv.gz",
    params:
        minRNACounts=lambda wc: config[wc.project]["configs"][wc.config]["minRNACounts"],
        minDNACounts=lambda wc: config[wc.project]["configs"][wc.config]["minDNACounts"],
    shell:
        """
        python {SCRIPTS_DIR}/count/merge_label.py --counts {input.counts} \
        --minRNACounts {params.minRNACounts} --minDNACounts {params.minDNACounts} \
        --assignment {input.association} \
        --output {output.counts} \
        --statistic {output.stats}
        """


rule make_master_tables:
    conda:
        "../envs/mpraflow_r.yaml"
    input:
        counts=lambda wc: expand(
            "results/{{project}}/assigned_counts/{{assignment}}/{{config}}/{{condition}}_{replicate}_merged_assigned_counts_{{sampling}}.tsv.gz",
            replicate=getReplicatesOfCondition(wc.project, wc.condition),
        ),
    output:
        statistic="results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_{sampling}_average_allreps_merged.tsv.gz",
        all="results/{project}/assigned_counts/{assignment}/{config}/{condition}_{sampling}_allreps_merged.tsv.gz",
        thresh="results/{project}/assigned_counts/{assignment}/{config}/{condition}_{sampling}_allreps_minThreshold_merged.tsv.gz",
    params:
        cond="{condition}",
        files=lambda wc: ",".join(
            expand(
                "results/{project}/assigned_counts/{assignment}/{config}/{condition}_{replicate}_merged_assigned_counts_{sampling}.tsv.gz",
                replicate=getReplicatesOfCondition(wc.project, wc.condition),
                project=wc.project,
                condition=wc.condition,
                assignment=wc.assignment,
                sampling=wc.sampling,
                config=wc.config,
            )
        ),
        replicates=lambda wc: ",".join(
            getReplicatesOfCondition(wc.project, wc.condition)
        ),
        thresh=lambda wc: config[wc.project]["configs"][wc.config]["bc_threshold"],
    shell:
        """
        Rscript {SCRIPTS_DIR}/count/make_master_tables.R \
        --condition {params.cond} \
        --threshold {params.thresh} \
        --files {params.files} \
        --replicates {params.replicates} \
        --output {output.all} \
        --output-all {output.thresh} \
        --statistic {output.statistic}
        """
