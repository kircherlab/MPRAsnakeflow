###################################
## BC overlap bweteen replicates ##
###################################


rule experiment_statistic_bc_overlap_run:
    """
    Get overlap of counts and barcodes between replicates.
    """
    conda:
        getCondaEnv("r.yaml")
    input:
        files=lambda wc: expand(
            getFinalCounts(
                wc.project, wc.config, wc.condition, wc.type, wc.raw_or_assigned
            ),
            project=wc.project,
            condition=wc.condition,
            config=wc.config,
            replicate=getReplicatesOfCondition(wc.project, wc.condition),
        ),
        script=getScript("count/BCCounts_betweenReplicates.R"),
    output:
        temp(
            "results/experiments/{project}/statistic/bc_overlap/{raw_or_assigned}/overlapBCandCounts.{condition}_{type}.{config}.tsv"
        ),
    params:
        input=lambda wc: ",".join(
            expand(
                getFinalCounts(
                    wc.project, wc.config, wc.condition, wc.type, wc.raw_or_assigned
                ),
                project=wc.project,
                condition=wc.condition,
                config=wc.config,
                replicate=getReplicatesOfCondition(wc.project, wc.condition),
            )
        ),
        cond="{condition}_{type}",
        replicates=lambda wc: ",".join(
            getReplicatesOfConditionType(wc.project, wc.condition, wc.type)
        ),
    log:
        temp(
            "results/logs/experiment/statistic/bc_overlap/run.{project}.{condition}.{type}.{config}.{raw_or_assigned}.log"
        ),
    shell:
        """
        Rscript {input.script} \
        --outfile {output} \
        --condition {params.cond} \
        --files {params.input} --replicates {params.replicates} &> {log}
        """


rule experiment_statistic_bc_overlap_combine_counts:
    """
    Combine overlap BC and count statistic into one file (raw counts).
    """
    conda:
        getCondaEnv("default.yaml")
    input:
        statistic=lambda wc: expand(
            "results/experiments/{{project}}/statistic/bc_overlap/counts/overlapBCandCounts.{condition}_{type}.{config}.tsv",
            type=["DNA", "RNA"],
            condition=getConditions(wc.project),
            config=wc.config,
        ),
    output:
        report(
            "results/experiments/{project}/statistic/bc_overlap.counts.{config}.tsv",
            caption="../../../report/bc_overlap.rst",
            category="{project}",
            subcategory="Barcodes",
            labels={
                "Analysis": "Overlap between replicates (not assigned)",
                "Configuration": "{config}",
                "Assignment": "-",
            },
        ),
    log:
        temp(
            "results/logs/experiment/statistic/bc_overlap/combine_counts.{project}.{config}.log"
        ),
    shell:
        """
        set +o pipefail;
        (
            cat {input.statistic[0]} | head -n 1;
            for i in {input.statistic}; do
                cat $i | tail -n +2
            done;
        ) > {output} 2> {log}
        """


rule experiment_statistic_bc_overlap_combine_assigned_counts:
    """
    Combine overlap BC and count statistic into one file (assigned counts).
    """
    conda:
        getCondaEnv("default.yaml")
    input:
        statistic=lambda wc: expand(
            "results/experiments/{{project}}/statistic/bc_overlap/assigned_counts/{{assignment}}/overlapBCandCounts.{condition}_{type}.{{config}}.tsv",
            type=["DNA", "RNA"],
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/experiments/{project}/statistic/bc_overlap.assigned_counts.{config}.{assignment}.tsv",
            caption="../../../report/bc_overlap_assignment.rst",
            category="{project}",
            subcategory="Barcodes",
            labels={
                "Analysis": "Overlap between replicates (assigned)",
                "Configuration": "{config}",
                "Assignment": "{assignment}",
            },
        ),
    log:
        temp(
            "results/logs/experiment/statistic/bc_overlap/combine_assigned_counts.{project}.{config}.{assignment}.log"
        ),
    shell:
        """
        set +o pipefail;
        (
            cat {input.statistic[0]} | head -n 1;
            for i in {input.statistic}; do
                cat $i | tail -n +2
            done;
        ) > {output} 2> {log}
        """
