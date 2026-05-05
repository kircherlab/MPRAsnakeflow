################################
### Correlation of BC Counts ###
################################


rule experiment_statistic_correlation_bc_counts:
    """
Calculate the correlation of the raw counts for each condition across replicates.
"""
    input:
        files=lambda wc: getMergedCounts(wc.project, wc.raw_or_assigned, wc.condition, wc.config)[0],
        script=getScript("count/plot_perBCCounts_correlation.R"),
    output:
        report(
            "results/experiments/{project}/statistic/barcode/{raw_or_assigned}/{condition}.{config}.barcode.DNA.pairwise.png",
            category="{project}",
            subcategory="BC correlation plots",
            labels={
                "Step": "{raw_or_assigned}",
                "Condition": "{condition}",
                "Configuration": "{config}",
                "Plot": "DNA",
            },
        ),
        report(
            "results/experiments/{project}/statistic/barcode/{raw_or_assigned}/{condition}.{config}.barcode.RNA.pairwise.png",
            category="{project}",
            subcategory="BC correlation plots",
            labels={
                "Step": "{raw_or_assigned}",
                "Condition": "{condition}",
                "Configuration": "{config}",
                "Plot": "RNA",
            },
        ),
        report(
            "results/experiments/{project}/statistic/barcode/{raw_or_assigned}/{condition}.{config}.barcode.Ratio.pairwise.png",
            category="{project}",
            subcategory="BC correlation plots",
            labels={
                "Step": "{raw_or_assigned}",
                "Condition": "{condition}",
                "Configuration": "{config}",
                "Plot": "Ratio",
            },
        ),
        temp("results/experiments/{project}/statistic/barcode/{raw_or_assigned}/{condition}.{config}.barcode.correlation.tsv"),
    log:
        temp(
            "results/logs/experiment/statistic/correlation/correlate_bc_counts.{project}.{condition}.{config}.{raw_or_assigned}.log"
        ),
    conda:
        getCondaEnv("r.yaml")
    params:
        replicates=lambda wc: ",".join(getMergedCounts(wc.project, wc.raw_or_assigned, wc.condition, wc.config)[1]),
        cond="{condition}",
        outdir="results/experiments/{project}/statistic/barcode/{raw_or_assigned}/{condition}.{config}",
        input=lambda wc: ",".join(getMergedCounts(wc.project, wc.raw_or_assigned, wc.condition, wc.config)[0]),
        minRNACounts=lambda wc: config["experiments"][wc.project]["configs"][wc.config]["filter"]["min_rna_counts"],
        minDNACounts=lambda wc: config["experiments"][wc.project]["configs"][wc.config]["filter"]["min_dna_counts"],
    shell:
        """
        Rscript {input.script} \
        --outdir {params.outdir} \
        --condition {params.cond} \
        --mindnacounts {params.minDNACounts} --minrnacounts {params.minRNACounts} \
        --files {params.input} --replicates {params.replicates} &> {log}
        """


rule experiment_statistic_correlation_bc_counts_hist:
    """
Generate histogram and boxplots of the raw counts for each condition across replicates.
"""
    input:
        files=lambda wc: getMergedCounts(wc.project, wc.raw_or_assigned, wc.condition, wc.config)[0],
        script=getScript("count/plot_perBCCounts_stats.R"),
    output:
        "results/experiments/{project}/statistic/barcode/{raw_or_assigned}/{condition}.{config}.DNA.perBarcode.png",
        "results/experiments/{project}/statistic/barcode/{raw_or_assigned}/{condition}.{config}.RNA.perBarcode.png",
    log:
        temp(
            "results/logs/experiment/statistic/correlation/correlate_bc_counts_hist.{project}.{condition}.{config}.{raw_or_assigned}.log"
        ),
    conda:
        getCondaEnv("r.yaml")
    params:
        replicates=lambda wc: ",".join(getMergedCounts(wc.project, wc.raw_or_assigned, wc.condition, wc.config)[1]),
        cond="{condition}",
        outdir="results/experiments/{project}/statistic/barcode/{raw_or_assigned}/{condition}.{config}",
        input=lambda wc: ",".join(getMergedCounts(wc.project, wc.raw_or_assigned, wc.condition, wc.config)[0]),
        minRNACounts=lambda wc: config["experiments"][wc.project]["configs"][wc.config]["filter"]["min_rna_counts"],
        minDNACounts=lambda wc: config["experiments"][wc.project]["configs"][wc.config]["filter"]["min_dna_counts"],
    shell:
        """
        Rscript {input.script} \
        --outdir {params.outdir} \
        --condition {params.cond} \
        --mindnacounts {params.minDNACounts} --minrnacounts {params.minRNACounts} \
        --files {params.input} --replicates {params.replicates} &> {log}
        """


rule experiment_statistic_correlation_combine_bc_raw:
    """
Combine the correlation of the raw counts for each condition across replicates into one table.
"""
    input:
        files=lambda wc: expand(
            "results/experiments/{{project}}/statistic/barcode/counts/{condition}.{{config}}.barcode.correlation.tsv",
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/experiments/{project}/statistic/statistic_bc_correlation_merged.{config}.tsv",
            caption="../../../report/bc_correlation.rst",
            category="{project}",
            subcategory="Barcode correlation",
            labels={
                "Configuration": "{config}",
                "Assignment": "-",
            },
        ),
    log:
        "results/logs/experiment/statistic/correlation/combine_bc_raw.{project}.{config}.log",
    conda:
        getCondaEnv("default.yaml")
    shell:
        """
        (
        cat {input.files[0]} | head -n 1;
        for i in {input.files}; do
            cat $i | tail -n +2;
        done;
        ) > {output} 2> {log}
        """


rule experiment_statistic_correlation_combine_bc_assigned:
    """
Combine the correlation of the assigned counts for each condition across replicates into one table.
"""
    input:
        files=lambda wc: expand(
            "results/experiments/{{project}}/statistic/barcode/assigned_counts/{{assignment}}/{condition}.{{config}}.barcode.correlation.tsv",
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/experiments/{project}/statistic/statistic_assigned_bc_correlation_merged.{assignment}.{config}.tsv",
            caption="../../../report/bc_correlation_assigned.rst",
            category="{project}",
            subcategory="Barcode correlation",
            labels={
                "Configuration": "{config}",
                "Assignment": "{assignment}",
            },
        ),
    log:
        temp("results/logs/experiment/statistic/correlation/combine_bc_assigned.{project}.{assignment}.{config}.log"),
    conda:
        getCondaEnv("default.yaml")
    shell:
        """
        (
        cat {input.files[0]} | head -n 1;
        for i in {input.files}; do
            cat $i | tail -n +2;
        done;
        ) > {output} 2> {log}
        """


#############################
### Correlation of Oligos ###
#############################


rule experiment_statistic_correlation_calculate:
    """
Calculate the correlation of oligos for each condition across replicates.
"""
    input:
        counts=lambda wc: expand(
            "results/experiments/{{project}}/assigned_counts/{{assignment}}/{{config}}/{{condition}}.{replicate}.merged_assigned_counts.tsv.gz",
            replicate=getReplicatesOfCondition(wc.project, wc.condition),
        ),
        label=(
            lambda wc: (
                config["experiments"][wc.project]["label_file"] if "label_file" in config["experiments"][wc.project] else []
            )
        ),
        script=getScript("count/plot_perInsertCounts_correlation.R"),
    output:
        report(
            "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}.DNA.pairwise.png",
            category="{project}",
            subcategory="Oligo correlation plots",
            labels={
                "Assignment": "{assignment}",
                "Condition": "{condition}",
                "Configuration": "{config}",
                "Plot": "DNA",
                "Threshold": "All",
            },
        ),
        report(
            "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}.RNA.pairwise.png",
            category="{project}",
            subcategory="Oligo correlation plots",
            labels={
                "Assignment": "{assignment}",
                "Condition": "{condition}",
                "Configuration": "{config}",
                "Plot": "RNA",
                "Threshold": "All",
            },
        ),
        report(
            "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}.Ratio.pairwise.png",
            category="{project}",
            subcategory="Oligo correlation plots",
            labels={
                "Assignment": "{assignment}",
                "Condition": "{condition}",
                "Configuration": "{config}",
                "Plot": "Ratio",
                "Threshold": "All",
            },
        ),
        report(
            "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}.DNA.pairwise.minThreshold.png",
            category="{project}",
            subcategory="Oligo correlation plots",
            labels=lambda wc: {
                "Assignment": "{assignment}",
                "Condition": "{condition}",
                "Configuration": "{config}",
                "Plot": "DNA",
                "Threshold": str(config["experiments"][wc.project]["configs"][wc.config]["filter"]["bc_threshold"]),
            },
        ),
        report(
            "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}.RNA.pairwise.minThreshold.png",
            category="{project}",
            subcategory="Oligo correlation plots",
            labels=lambda wc: {
                "Assignment": "{assignment}",
                "Condition": "{condition}",
                "Configuration": "{config}",
                "Plot": "RNA",
                "Threshold": str(config["experiments"][wc.project]["configs"][wc.config]["filter"]["bc_threshold"]),
            },
        ),
        report(
            "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}.Ratio.pairwise.minThreshold.png",
            category="{project}",
            subcategory="Oligo correlation plots",
            labels=lambda wc: {
                "Assignment": "{assignment}",
                "Condition": "{condition}",
                "Configuration": "{config}",
                "Plot": "Ratio",
                "Threshold": str(config["experiments"][wc.project]["configs"][wc.config]["filter"]["bc_threshold"]),
            },
        ),
        temp("results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}.correlation.tsv"),
        temp(
            "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}.correlation.minThreshold.tsv"
        ),
    log:
        temp("results/logs/experiment/statistic/correlation/calculate.{project}.{condition}.{config}.{assignment}.log"),
    conda:
        getCondaEnv("r.yaml")
    params:
        cond="{condition}",
        files=lambda wc: ",".join(
            expand(
                "results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}.{replicate}.merged_assigned_counts.tsv.gz",
                replicate=getReplicatesOfCondition(wc.project, wc.condition),
                project=wc.project,
                condition=wc.condition,
                assignment=wc.assignment,
                config=wc.config,
            )
        ),
        replicates=lambda wc: ",".join(getReplicatesOfCondition(wc.project, wc.condition)),
        thresh=lambda wc: config["experiments"][wc.project]["configs"][wc.config]["filter"]["bc_threshold"],
        outdir="results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}",
        label=(
            lambda wc: (
                "--label %s" % config["experiments"][wc.project]["label_file"]
                if "label_file" in config["experiments"][wc.project]
                else ""
            )
        ),
    shell:
        """
        Rscript {input.script} \
        --condition {params.cond} \
        {params.label} \
        --files {params.files} \
        --replicates {params.replicates} \
        --threshold {params.thresh} \
        --outdir {params.outdir} &> {log}
        """


rule experiment_statistic_correlation_hist_box_plots:
    """
Generate histogram and boxplots of the oligos for each condition across replicates.
"""
    input:
        counts=lambda wc: expand(
            "results/experiments/{{project}}/assigned_counts/{{assignment}}/{{config}}/{{condition}}.{replicate}.merged_assigned_counts.tsv.gz",
            replicate=getReplicatesOfCondition(wc.project, wc.condition),
        ),
        label=(
            lambda wc: (
                config["experiments"][wc.project]["label_file"] if "label_file" in config["experiments"][wc.project] else []
            )
        ),
        script=getScript("count/plot_perInsertCounts_stats.R"),
    output:
        "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}.barcodesPerInsert.png",
        "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}.group_barcodesPerInsert_box.png",
        "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}.group_barcodesPerInsert_box_minThreshold.png",
        "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}.dna_vs_rna.png",
        "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}.dna_vs_rna_minThreshold.png",
    log:
        temp("results/logs/experiment/statistic/correlation/hist_box_plots.{project}.{condition}.{config}.{assignment}.log"),
    conda:
        getCondaEnv("r.yaml")
    params:
        cond="{condition}",
        files=lambda wc: ",".join(
            expand(
                "results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}.{replicate}.merged_assigned_counts.tsv.gz",
                replicate=getReplicatesOfCondition(wc.project, wc.condition),
                project=wc.project,
                condition=wc.condition,
                assignment=wc.assignment,
                config=wc.config,
            )
        ),
        replicates=lambda wc: ",".join(getReplicatesOfCondition(wc.project, wc.condition)),
        thresh=lambda wc: config["experiments"][wc.project]["configs"][wc.config]["filter"]["bc_threshold"],
        outdir="results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}",
        label=(
            lambda wc: (
                "--label %s" % config["experiments"][wc.project]["label_file"]
                if "label_file" in config["experiments"][wc.project]
                else ""
            )
        ),
    shell:
        """
        Rscript {input.script} \
        --condition {params.cond} \
        {params.label} \
        --files {params.files} \
        --replicates {params.replicates} \
        --threshold {params.thresh} \
        --outdir {params.outdir} &> {log}
        """


rule experiment_statistic_correlation_combine_oligo:
    """
Combine the correlation of oligos for each condition across replicates into one table.
"""
    input:
        correlation=lambda wc: expand(
            "results/experiments/{{project}}/statistic/assigned_counts/{{assignment}}/{{config}}/{condition}.correlation.tsv",
            condition=getConditions(wc.project),
        ),
        correlation_thresh=lambda wc: expand(
            "results/experiments/{{project}}/statistic/assigned_counts/{{assignment}}/{{config}}/{condition}.correlation.minThreshold.tsv",
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/experiments/{project}/statistic/statistic_oligo_correlation_merged.{assignment}.{config}.tsv",
            caption="../../../report/oligo_correlation.rst",
            category="{project}",
            subcategory="Oligo correlation",
            labels={
                "Analysis": "Oligo correlation",
                "Configuration": "{config}",
                "Assignment": "{assignment}",
            },
        ),
    log:
        temp("results/logs/experiment/statistic/correlation/combine_oligo.{project}.{config}.{assignment}.log"),
    conda:
        getCondaEnv("default.yaml")
    params:
        thresh=lambda wc: config["experiments"][wc.project]["configs"][wc.config]["filter"]["bc_threshold"],
    shell:
        """
        set +o pipefail;
        (
        cat {input.correlation[0]} | head -n 1 | awk -v 'OFS=\\t' '{{print $0,"threshold (min {params.thresh})"}}';
        for i in {input.correlation}; do
            cat $i | tail -n +2 | awk -v 'OFS=\\t' '{{print $0,"False"}}'
        done;
        for i in {input.correlation_thresh}; do
            cat $i | tail -n +2 | awk -v 'OFS=\\t' '{{print $0,"True"}}'
        done;
        ) > {output} 2> {log}
        """
