# --use-pseudo-counts
rule variants_generateVariantTable:
    input:
        variant_definition=lambda wc: getVariants(wc.project)["map"],
        counts="results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}.{replicate}.merged_assigned_counts.tsv.gz",
        script=getScript("variants/generateVariantTable.py"),
    output:
        "results/experiments/{project}/variants/{assignment}/{config}/{condition}.{replicate}.variantTable.tsv.gz",
    log:
        temp(
            "results/logs/experiments/variants/generateVariantTable.{project}.{assignment}.{config}.{condition}.{replicate}.log"
        ),
    conda:
        getCondaEnv("python3.yaml")
    shell:
        """
        python {input.script} \
        --counts {input.counts} \
        --declaration {input.variant_definition} \
        --output {output} &> {log}
        """


rule variants_MasterTable:
    input:
        variants=lambda wc: expand(
            "results/experiments/{{project}}/variants/{{assignment}}/{{config}}/{{condition}}.{replicate}.variantTable.tsv.gz",
            replicate=getReplicatesOfCondition(wc.project, wc.condition),
        ),
        script=getScript("variants/generateMasterVariantTable.py"),
    output:
        "results/experiments/{project}/variants/{assignment}/{config}/{condition}.variantTable.tsv.gz",
    log:
        temp("results/logs/experiments/variants/MasterTable.{project}.{assignment}/{config}.{condition}.log"),
    conda:
        getCondaEnv("python3.yaml")
    params:
        input=lambda wc: " ".join(
            [
                "--input %s" % i
                for i in expand(
                    "results/experiments/{project}/variants/{assignment}/{config}/{condition}.{replicate}.variantTable.tsv.gz",
                    replicate=getReplicatesOfCondition(wc.project, wc.condition),
                    project=wc.project,
                    assignment=wc.assignment,
                    config=wc.config,
                    condition=wc.condition,
                )
            ]
        ),
        minRNACounts=lambda wc: config["experiments"][wc.project]["configs"][wc.config]["filter"]["min_rna_counts"],
        minDNACounts=lambda wc: config["experiments"][wc.project]["configs"][wc.config]["filter"]["min_dna_counts"],
    shell:
        """
        python {input.script} \
        {params.input} \
        --minRNACounts {params.minRNACounts} \
        --minDNACounts {params.minDNACounts} \
        --output {output} &> {log}
        """


rule variants_correlate:
    input:
        counts=lambda wc: expand(
            "results/experiments/{{project}}/variants/{{assignment}}/{{config}}/{{condition}}.{replicate}.variantTable.tsv.gz",
            replicate=getReplicatesOfCondition(wc.project, wc.condition),
        ),
        script=getScript("variants/correlateVariantTables.py"),
    output:
        "results/experiments/{project}/statistic/variants/{assignment}/{config}/{condition}/{condition}.correlation_variantTable_minBC{threshold}.tsv.gz",
    log:
        temp(
            "results/logs/experiments/variants/correlate.{project}.{assignment}.{config}.{condition}.{condition}.{threshold}.log"
        ),
    conda:
        getCondaEnv("python3.yaml")
    params:
        cond="{condition}",
        threshold="{threshold}",
        tables=lambda wc: " ".join(
            expand(
                "--variants {replicate} results/experiments/{project}/variants/{assignment}/{config}/{condition}.{replicate}.variantTable.tsv.gz",
                project=wc.project,
                assignment=wc.assignment,
                config=wc.config,
                condition=wc.condition,
                replicate=getReplicatesOfCondition(wc.project, wc.condition),
            )
        ),
    shell:
        """
        python {input.script} \
        --condition {params.cond} \
        {params.tables} \
        --bc-threshold {params.threshold} \
        --output {output} &> {log}
        """


rule variants_combineVariantCorrelationTables:
    input:
        correlation=lambda wc: expand(
            "results/experiments/{{project}}/statistic/variants/{{assignment}}/{{config}}/{condition}/{condition}.correlation_variantTable_minBC{threshold}.tsv.gz",
            condition=getConditions(wc.project),
            threshold=getVariantsBCThreshold(wc.project),
        ),
    output:
        report(
            "results/experiments/{project}/statistic/variants/{assignment}/{config}/correlation_variantTable.tsv",
            caption="../report/variants/correlation.rst",
            category="{project}",
            subcategory="Variants",
            labels={
                "Analysis": "Correlation",
                "Configuration": "{config}",
                "Assignment": "{assignment}",
            },
        ),
    log:
        temp("results/logs/experiments/variants/combineVariantCorrelationTables.{project}.{assignment}.{config}.log"),
    conda:
        getCondaEnv("default.yaml")
    shell:
        """
        set +o pipefail;
        (
        zcat {input.correlation[0]} | head -n 1;
        for i in {input.correlation}; do
            zcat $i | tail -n +2;
        done;
        ) > {output} 2> {log}
        """
