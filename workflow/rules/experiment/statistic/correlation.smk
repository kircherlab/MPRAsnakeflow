################################
### Correlation of BC Counts ###
################################


rule experiment_statistic_correlation_bc_counts:
    conda:
        getCondaEnv("r.yaml")
    input:
        files=lambda wc: getMergedCounts(
            wc.project, wc.raw_or_assigned, wc.condition, wc.config
        )[0],
        script=getScript("count/plot_perBCCounts_correlation.R"),
    output:
        report(
            "results/experiments/{project}/statistic/barcode/{raw_or_assigned}/{condition}_{config}_barcode_DNA_pairwise.png",
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
            "results/experiments/{project}/statistic/barcode/{raw_or_assigned}/{condition}_{config}_barcode_RNA_pairwise.png",
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
            "results/experiments/{project}/statistic/barcode/{raw_or_assigned}/{condition}_{config}_barcode_Ratio_pairwise.png",
            category="{project}",
            subcategory="BC correlation plots",
            labels={
                "Step": "{raw_or_assigned}",
                "Condition": "{condition}",
                "Configuration": "{config}",
                "Plot": "Ratio",
            },
        ),
        temp(
            "results/experiments/{project}/statistic/barcode/{raw_or_assigned}/{condition}_{config}_barcode_correlation.tsv"
        ),
    params:
        replicates=lambda wc: ",".join(
            getMergedCounts(wc.project, wc.raw_or_assigned, wc.condition, wc.config)[1]
        ),
        cond="{condition}",
        outdir="results/experiments/{project}/statistic/barcode/{raw_or_assigned}/{condition}_{config}",
        input=lambda wc: ",".join(
            getMergedCounts(wc.project, wc.raw_or_assigned, wc.condition, wc.config)[0]
        ),
        minRNACounts=lambda wc: config["experiments"][wc.project]["configs"][
            wc.config
        ]["filter"]["min_rna_counts"],
        minDNACounts=lambda wc: config["experiments"][wc.project]["configs"][
            wc.config
        ]["filter"]["min_dna_counts"],
    log:
        temp(
            "results/logs/experiment/statistic/correlation/correlate_bc_counts.{project}.{condition}.{config}.{raw_or_assigned}.log"
        ),
    shell:
        """
        Rscript {input.script} \
        --outdir {params.outdir} \
        --condition {params.cond} \
        --mindnacounts {params.minDNACounts} --minrnacounts {params.minRNACounts} \
        --files {params.input} --replicates {params.replicates} &> {log}
        """


rule experiment_statistic_correlation_bc_counts_hist:
    conda:
        getCondaEnv("r.yaml")
    input:
        files=lambda wc: getMergedCounts(
            wc.project, wc.raw_or_assigned, wc.condition, wc.config
        )[0],
        script=getScript("count/plot_perBCCounts_stats.R"),
    output:
        "results/experiments/{project}/statistic/barcode/{raw_or_assigned}/{condition}_{config}_DNA_perBarcode.png",
        "results/experiments/{project}/statistic/barcode/{raw_or_assigned}/{condition}_{config}_RNA_perBarcode.png",
    params:
        replicates=lambda wc: ",".join(
            getMergedCounts(wc.project, wc.raw_or_assigned, wc.condition, wc.config)[1]
        ),
        cond="{condition}",
        outdir="results/experiments/{project}/statistic/barcode/{raw_or_assigned}/{condition}_{config}",
        input=lambda wc: ",".join(
            getMergedCounts(wc.project, wc.raw_or_assigned, wc.condition, wc.config)[0]
        ),
        minRNACounts=lambda wc: config["experiments"][wc.project]["configs"][
            wc.config
        ]["filter"]["min_rna_counts"],
        minDNACounts=lambda wc: config["experiments"][wc.project]["configs"][
            wc.config
        ]["filter"]["min_dna_counts"],
    log:
        temp(
            "results/logs/experiment/statistic/correlation/correlate_bc_counts_hist.{project}.{condition}.{config}.{raw_or_assigned}.log"
        ),
    shell:
        """
        Rscript {input.script} \
        --outdir {params.outdir} \
        --condition {params.cond} \
        --mindnacounts {params.minDNACounts} --minrnacounts {params.minRNACounts} \
        --files {params.input} --replicates {params.replicates} &> {log}
        """


rule experiment_statistic_correlation_combine_bc_raw:
    conda:
        getCondaEnv("default.yaml")
    input:
        files=lambda wc: expand(
            "results/experiments/{{project}}/statistic/barcode/counts/{condition}_{{config}}_barcode_correlation.tsv",
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/experiments/{project}/statistic/statistic_bc_correlation_merged_{config}.tsv",
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
    conda:
        getCondaEnv("default.yaml")
    input:
        files=lambda wc: expand(
            "results/experiments/{{project}}/statistic/barcode/assigned_counts/{{assignment}}/{condition}_{{config}}_barcode_correlation.tsv",
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/experiments/{project}/statistic/statistic_assigned_bc_correlation_merged_{assignment}_{config}.tsv",
            caption="../../../report/bc_correlation_assigned.rst",
            category="{project}",
            subcategory="Barcode correlation",
            labels={
                "Configuration": "{config}",
                "Assignment": "{assignment}",
            },
        ),
    log:
        temp(
            "results/logs/experiment/statistic/correlation/combine_bc_assigned.{project}.{assignment}_{config}.log"
        ),
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
    conda:
        getCondaEnv("r.yaml")
    input:
        counts=lambda wc: expand(
            "results/experiments/{{project}}/assigned_counts/{{assignment}}/{{config}}/{{condition}}_{replicate}_merged_assigned_counts.tsv.gz",
            replicate=getReplicatesOfCondition(wc.project, wc.condition),
        ),
        label=(
            lambda wc: (
                config["experiments"][wc.project]["label_file"]
                if "label_file" in config["experiments"][wc.project]
                else []
            )
        ),
        script=getScript("count/plot_perInsertCounts_correlation.R"),
    output:
        report(
            "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_DNA_pairwise.png",
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
            "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_RNA_pairwise.png",
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
            "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_Ratio_pairwise.png",
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
            "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_DNA_pairwise_minThreshold.png",
            category="{project}",
            subcategory="Oligo correlation plots",
            labels=lambda wc: {
                "Assignment": "{assignment}",
                "Condition": "{condition}",
                "Configuration": "{config}",
                "Plot": "DNA",
                "Threshold": str(
                    config["experiments"][wc.project]["configs"][wc.config]["filter"][
                        "bc_threshold"
                    ]
                ),
            },
        ),
        report(
            "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_RNA_pairwise_minThreshold.png",
            category="{project}",
            subcategory="Oligo correlation plots",
            labels=lambda wc: {
                "Assignment": "{assignment}",
                "Condition": "{condition}",
                "Configuration": "{config}",
                "Plot": "RNA",
                "Threshold": str(
                    config["experiments"][wc.project]["configs"][wc.config]["filter"][
                        "bc_threshold"
                    ]
                ),
            },
        ),
        report(
            "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_Ratio_pairwise_minThreshold.png",
            category="{project}",
            subcategory="Oligo correlation plots",
            labels=lambda wc: {
                "Assignment": "{assignment}",
                "Condition": "{condition}",
                "Configuration": "{config}",
                "Plot": "Ratio",
                "Threshold": str(
                    config["experiments"][wc.project]["configs"][wc.config]["filter"][
                        "bc_threshold"
                    ]
                ),
            },
        ),
        temp(
            "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_correlation.tsv"
        ),
        temp(
            "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_correlation_minThreshold.tsv"
        ),
    params:
        cond="{condition}",
        files=lambda wc: ",".join(
            expand(
                "results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}_{replicate}_merged_assigned_counts.tsv.gz",
                replicate=getReplicatesOfCondition(wc.project, wc.condition),
                project=wc.project,
                condition=wc.condition,
                assignment=wc.assignment,
                config=wc.config,
            )
        ),
        replicates=lambda wc: ",".join(
            getReplicatesOfCondition(wc.project, wc.condition)
        ),
        thresh=lambda wc: config["experiments"][wc.project]["configs"][wc.config][
            "filter"
        ]["bc_threshold"],
        outdir="results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}",
        label=(
            lambda wc: (
                "--label %s" % config["experiments"][wc.project]["label_file"]
                if "label_file" in config["experiments"][wc.project]
                else ""
            )
        ),
    log:
        temp(
            "results/logs/experiment/statistic/correlation/calculate.{project}.{condition}.{config}.{assignment}.log"
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
    conda:
        getCondaEnv("r.yaml")
    input:
        counts=lambda wc: expand(
            "results/experiments/{{project}}/assigned_counts/{{assignment}}/{{config}}/{{condition}}_{replicate}_merged_assigned_counts.tsv.gz",
            replicate=getReplicatesOfCondition(wc.project, wc.condition),
        ),
        label=(
            lambda wc: (
                config["experiments"][wc.project]["label_file"]
                if "label_file" in config["experiments"][wc.project]
                else []
            )
        ),
        script=getScript("count/plot_perInsertCounts_stats.R"),
    output:
        "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_barcodesPerInsert.png",
        "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_group_barcodesPerInsert_box.png",
        "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_group_barcodesPerInsert_box_minThreshold.png",
        "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_dna_vs_rna.png",
        "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_dna_vs_rna_minThreshold.png",
    params:
        cond="{condition}",
        files=lambda wc: ",".join(
            expand(
                "results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}_{replicate}_merged_assigned_counts.tsv.gz",
                replicate=getReplicatesOfCondition(wc.project, wc.condition),
                project=wc.project,
                condition=wc.condition,
                assignment=wc.assignment,
                config=wc.config,
            )
        ),
        replicates=lambda wc: ",".join(
            getReplicatesOfCondition(wc.project, wc.condition)
        ),
        thresh=lambda wc: config["experiments"][wc.project]["configs"][wc.config][
            "filter"
        ]["bc_threshold"],
        outdir="results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}",
        label=(
            lambda wc: (
                "--label %s" % config["experiments"][wc.project]["label_file"]
                if "label_file" in config["experiments"][wc.project]
                else ""
            )
        ),
    log:
        temp(
            "results/logs/experiment/statistic/correlation/hist_box_plots.{project}.{condition}.{config}.{assignment}.log"
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
    conda:
        getCondaEnv("default.yaml")
    input:
        correlation=lambda wc: expand(
            "results/experiments/{{project}}/statistic/assigned_counts/{{assignment}}/{{config}}/{condition}_correlation.tsv",
            condition=getConditions(wc.project),
        ),
        correlation_thresh=lambda wc: expand(
            "results/experiments/{{project}}/statistic/assigned_counts/{{assignment}}/{{config}}/{condition}_correlation_minThreshold.tsv",
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/experiments/{project}/statistic/statistic_oligo_correlation_merged_{assignment}_{config}.tsv",
            caption="../../../report/oligo_correlation.rst",
            category="{project}",
            subcategory="Oligo correlation",
            labels={
                "Analysis": "Oligo correlation",
                "Configuration": "{config}",
                "Assignment": "{assignment}",
            },
        ),
    params:
        thresh=lambda wc: config["experiments"][wc.project]["configs"][wc.config][
            "filter"
        ]["bc_threshold"],
    log:
        temp(
            "results/logs/experiment/statistic/correlation/combine_oligo.{project}.{config}.{assignment}.log"
        ),
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
