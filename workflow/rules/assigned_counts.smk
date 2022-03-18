##################################
### Assign BCs and afterwards ###
##################################


rule assignBarcodes:
    """
    Assign RNA and DNA barcodes seperately to make the statistic for assigned
    """
    conda:
        "../envs/python3.yaml"
    input:
        counts="results/experiments/{project}/counts/{condition}_{replicate}_{type}_final_counts_{sampling}.tsv.gz",
        association=lambda wc: getAssignmentFile(wc.project, wc.assignment),
    output:
        counts="results/experiments/{project}/assigned_counts/{assignment}/{condition}_{replicate}_{type}_final_counts_{sampling}.tsv.gz",
        stats="results/experiments/{project}/stats/assigned_counts/{assignment}/{condition}_{replicate}_{type}_{sampling}.statistic.tsv.gz",
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
        "envs/python3.yaml"
    input:
        lambda wc: getAssignmentFile(wc.project, wc.assignment),
    output:
        "results/experiments/{project}/assigned_counts/{assignment}/assignment.pickle",
    shell:
        """
        python {SCRIPTS_DIR}/count/create_pickle.py -i {input} -o {output}
        """


rule dna_rna_merge:
    conda:
        "../envs/python3.yaml"
    input:
        counts=lambda wc: expand(
            "results/experiments/{{project}}/counts/merged/{mergeType}/{{condition}}_{{replicate}}_merged_counts_{{sampling}}.tsv.gz",
        mergeType="withoutZeros"
            if config["experiments"][wc.project]["configs"][wc.config]["minRNACounts"]
            > 0
            and config["experiments"][wc.project]["configs"][wc.config]["minDNACounts"]
            > 0
            else "withZeros",
        ),
        association=lambda wc: getAssignmentFile(wc.project, wc.assignment),
    output:
        counts="results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}_{replicate}_merged_assigned_counts_{sampling}.tsv.gz",
        stats="results/experiments/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_{replicate}_merged_assigned_counts_{sampling}.statistic.tsv.gz",
    params:
        minRNACounts=lambda wc: config["experiments"][wc.project]["configs"][
            wc.config
        ]["minRNACounts"],
        minDNACounts=lambda wc: config["experiments"][wc.project]["configs"][
            wc.config
        ]["minDNACounts"],
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
        "../envs/r.yaml"
    input:
        counts=lambda wc: expand(
            "results/experiments/{{project}}/assigned_counts/{{assignment}}/{{config}}/{{condition}}_{replicate}_merged_assigned_counts_{{sampling}}.tsv.gz",
            replicate=getReplicatesOfCondition(wc.project, wc.condition),
        ),
    output:
        statistic="results/experiments/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_{sampling}_average_allreps_merged.tsv.gz",
        all="results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}_{sampling}_allreps_merged.tsv.gz",
        thresh="results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}_{sampling}_allreps_minThreshold_merged.tsv.gz",
    params:
        cond="{condition}",
        files=lambda wc: ",".join(
            expand(
                "results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}_{replicate}_merged_assigned_counts_{sampling}.tsv.gz",
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
        thresh=lambda wc: config["experiments"][wc.project]["configs"][wc.config][
            "bc_threshold"
        ],
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
