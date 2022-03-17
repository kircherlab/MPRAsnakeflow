###################################
## BC overlap bweteen replicates ##
###################################


rule statistic_overlapBCs:
    conda:
        "../../envs/mpraflow_r.yaml"
    input:
        lambda wc: expand(
            "results/experiments/{{project}}/{{raw_or_assigned}}/{{condition}}_{replicate}_{{type}}_final_counts_{{sampling}}.tsv.gz",
            replicate=getReplicatesOfCondition(wc.project, wc.condition),
        ),
    output:
        "results/experiments/{project}/stats/{raw_or_assigned}/overlapBCandCounts_{condition}_{type}_{sampling}.tsv",
    params:
        input=lambda wc: ",".join(
            expand(
                "results/experiments/{project}/{raw_or_assigned}/{condition}_{replicate}_{type}_final_counts_{sampling}.tsv.gz",
                project=wc.project,
                condition=wc.condition,
                raw_or_assigned=wc.raw_or_assigned,
                type=wc.type,
                sampling=wc.sampling,
                replicate=getReplicatesOfCondition(wc.project, wc.condition),
            )
        ),
        cond="{condition}_{type}",
        replicates=lambda wc: ",".join(
            getReplicatesOfCondition(wc.project, wc.condition)
        ),
    shell:
        """
        Rscript {SCRIPTS_DIR}/count/BCCounts_betweenReplicates.R \
        --outfile {output} \
        --condition {params.cond} \
        --files {params.input} --replicates {params.replicates}
        """


rule statistic_combine_overlapBCs_stats_raw:
    input:
        stats=lambda wc: expand(
            "results/experiments/{{project}}/stats/counts/overlapBCandCounts_{condition}_{type}_{sampling}.tsv",
            type=["DNA", "RNA"],
            condition=getConditions(wc.project),
            sampling=wc.sampling,
        ),
    output:
        report(
            "results/experiments/{project}/stats/statistic_overlapBCs_counts_{sampling}.tsv",
            caption="../../report/bc_overlap.rst",
            category="{project}",
            subcategory="Barcodes",
        ),
    shell:
        """
        set +o pipefail;
        (
            cat {input.stats[0]} | head -n 1;
            for i in {input.stats}; do
                cat $i | tail -n +2
            done;
        ) > {output}
        """


rule statistic_combine_overlapBCs_stats_assigned:
    input:
        stats=lambda wc: expand(
            "results/experiments/{{project}}/stats/assigned_counts/{{assignment}}/overlapBCandCounts_{condition}_{type}_{{sampling}}.tsv",
            type=["DNA", "RNA"],
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/experiments/{project}/stats/statistic_overlapBCs_assigned_counts_{assignment}_{sampling}.tsv",
            caption="../../report/bc_overlap_assignment.rst",
            category="{project}",
            subcategory="Barcodes",
        ),
    shell:
        """
        set +o pipefail;
        (
            cat {input.stats[0]} | head -n 1;
            for i in {input.stats}; do
                cat $i | tail -n +2
            done;
        ) > {output}
        """
