

rule generateVariantTable:
    conda:
        "../envs/mpraflow_py36.yaml"
    input:
        variant_definition=lambda wc: getVariants(wc.project)["map"],
        counts="results/{project}/assigned_counts/{assignment}/{config}/{condition}_{replicate}_merged_assigned_counts.tsv.gz",
    output:
        "results/{project}/variants/{assignment}/{config}/{condition}_{replicate}_variantTable.tsv.gz",
    shell:
        """
        python {SCRIPTS_DIR}/variants/generateVariantTable.py \
        --counts {input.counts} \
        --declaration {input.variant_definition} \
        --output {output}
        """


rule correlate_variants:
    conda:
        "../envs/mpraflow_py36.yaml"
    input:
        counts=lambda wc: expand(
            "results/{{project}}/variants/{{assignment}}/{{config}}/{{condition}}_{replicate}_variantTable.tsv.gz",
            replicate=getReplicatesOfCondition(wc.project, wc.condition),
        ),
    output:
        "results/{project}/stats/variants/{assignment}/{config}/{condition}/{condition}_correlation_variantTable_minBC{threshold}.tsv.gz",
    params:
        cond="{condition}",
        threshold="{threshold}",
        tables=lambda wc: " ".join(
            expand(
                "--variants {replicate} results/{project}/variants/{assignment}/{config}/{condition}_{replicate}_variantTable.tsv.gz",
                project=wc.project,
                assignment=wc.assignment,
                config=wc.config,
                condition=wc.condition,
                replicate=getReplicatesOfCondition(wc.project, wc.condition),
            )
        ),
    shell:
        """
        python {SCRIPTS_DIR}/variants/correlateVariantTables.py \
        --condition {params.cond} \
        {params.tables} \
        --bc-threshold {params.threshold} \
        --output {output}
        """


rule combineVariantCorrelationTables:
    input:
        correlation=lambda wc: expand(
            "results/{{project}}/stats/variants/{{assignment}}/{{config}}/{condition}/{condition}_correlation_variantTable_minBC{threshold}.tsv.gz",
            condition=getConditions(wc.project),
            threshold=getVariantsBCThreshold(wc.project),
        ),
    output:
        report(
            "results/{project}/stats/variants/{assignment}/{config}/correlation_variantTable.tsv",
            caption="../report/variants/correlation.rst",
            category="{project}",
            subcategory="Variants",
        ),

    shell:
        """
        set +o pipefail;
        (
        zcat {input.correlation[0]} | head -n 1;
        for i in {input.correlation}; do
            zcat $i | tail -n +2;
        done;
        ) > {output}
        """
