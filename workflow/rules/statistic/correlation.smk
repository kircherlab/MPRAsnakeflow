################################
### Correlation of BC Counts ###
################################


rule statistic_correlate_BC_counts:
    conda:
        "../../envs/r.yaml"
    input:
        files=lambda wc: getMergedCounts(
            wc.project, wc.raw_or_assigned, wc.condition, wc.config
        )[0],
        script=getScript("count/plot_perBCCounts_correlation.R"),
    output:
        "results/experiments/{project}/stats/barcode/{raw_or_assigned}/{condition}_{config}_barcode_DNA_pairwise.png",
        "results/experiments/{project}/stats/barcode/{raw_or_assigned}/{condition}_{config}_barcode_RNA_pairwise.png",
        "results/experiments/{project}/stats/barcode/{raw_or_assigned}/{condition}_{config}_barcode_Ratio_pairwise.png",
        "results/experiments/{project}/stats/barcode/{raw_or_assigned}/{condition}_{config}_barcode_correlation.tsv",
        "results/experiments/{project}/stats/barcode/{raw_or_assigned}/{condition}_{config}_DNA_perBarcode.png",
        "results/experiments/{project}/stats/barcode/{raw_or_assigned}/{condition}_{config}_RNA_perBarcode.png",
    params:
        replicates=lambda wc: ",".join(
            getMergedCounts(wc.project, wc.raw_or_assigned, wc.condition, wc.config)[1]
        ),
        cond="{condition}",
        outdir="results/experiments/{project}/stats/barcode/{raw_or_assigned}/{condition}_{config}",
        input=lambda wc: ",".join(
            getMergedCounts(wc.project, wc.raw_or_assigned, wc.condition, wc.config)[0]
        ),
        minRNACounts=lambda wc: config["experiments"][wc.project]["configs"][
            wc.config
        ]["filter"]["RNA"]["minCounts"],
        minDNACounts=lambda wc: config["experiments"][wc.project]["configs"][
            wc.config
        ]["filter"]["DNA"]["minCounts"],
    log:
        "logs/experiments/{project}/stats/barcode/{raw_or_assigned}/statistic_correlate_BC_counts.{condition}_{config}.log",
    shell:
        """
        Rscript {input.script} \
        --outdir {params.outdir} \
        --condition {params.cond} \
        --mindnacounts {params.minDNACounts} --minrnacounts {params.minRNACounts} \
        --files {params.input} --replicates {params.replicates} &> {log}
        """


rule statistic_combine_bc_correlation_raw:
    conda:
        "../../envs/default.yaml"
    input:
        files=lambda wc: expand(
            "results/experiments/{{project}}/stats/barcode/counts/{condition}_{{config}}_barcode_correlation.tsv",
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/experiments/{project}/stats/statistic_bc_correlation_merged_{config}.tsv",
            caption="../../report/bc_correlation.rst",
            category="{project}",
            subcategory="Barcodes",
        ),
    log:
        "logs/experiments/{project}/stats/statistic_combine_bc_correlation_raw.{config}.log",
    shell:
        """
        (
        cat {input.files[0]} | head -n 1;
        for i in {input.files}; do
            cat $i | tail -n +2;
        done;
        ) > {output}
        """


rule statistic_combine_bc_correlation_assigned:
    conda:
        "../../envs/default.yaml"
    input:
        files=lambda wc: expand(
            "results/experiments/{{project}}/stats/barcode/assigned_counts/{{assignment}}/{condition}_{{config}}_barcode_correlation.tsv",
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/experiments/{project}/stats/statistic_assigned_bc_correlation_merged_{assignment}_{config}.tsv",
            caption="../../report/bc_correlation_assigned.rst",
            category="{project}",
            subcategory="Barcodes",
        ),
    log:
        "logs/experiments/{project}/stats/statistic_combine_bc_correlation_assigned{assignment}_{config}.log",
    shell:
        """
        (
        cat {input.files[0]} | head -n 1;
        for i in {input.files}; do
            cat $i | tail -n +2;
        done;
        ) > {output}
        """


#############################
### Correlation of Oligos ###
#############################


rule statistic_calc_correlations:
    conda:
        "../../envs/r.yaml"
    input:
        counts=lambda wc: expand(
            "results/experiments/{{project}}/assigned_counts/{{assignment}}/{{config}}/{{condition}}_{replicate}_merged_assigned_counts.tsv.gz",
            replicate=getReplicatesOfCondition(wc.project, wc.condition),
        ),
        label=(
            lambda wc: config["experiments"][wc.project]["label_file"]
            if "label_file" in config["experiments"][wc.project]
            else []
        ),
        script=getScript("count/plot_perInsertCounts_correlation.R"),
    output:
        "results/experiments/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_all_barcodesPerInsert_box.png",
        "results/experiments/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_all_barcodesPerInsert_box_minThreshold.png",
        "results/experiments/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_DNA_pairwise.png",
        "results/experiments/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_DNA_pairwise_minThreshold.png",
        "results/experiments/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_group_barcodesPerInsert_box.png",
        "results/experiments/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_group_barcodesPerInsert_box_minThreshold.png",
        "results/experiments/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_Ratio_pairwise.png",
        "results/experiments/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_Ratio_pairwise_minThreshold.png",
        "results/experiments/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_RNA_pairwise.png",
        "results/experiments/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_RNA_pairwise_minThreshold.png",
        "results/experiments/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_correlation.tsv",
        "results/experiments/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_correlation_minThreshold.tsv",
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
            "bc_threshold"
        ],
        outdir="results/experiments/{project}/stats/assigned_counts/{assignment}/{config}/{condition}",
        label=(
            lambda wc: "--label %s" % config["experiments"][wc.project]["label_file"]
            if "label_file" in config["experiments"][wc.project]
            else ""
        ),
    log:
        "logs/experiments/{project}/stats/assigned_counts/{assignment}/{config}/statistic_calc_correlations.{condition}.log",
    shell:
        """
        Rscript {input.script} \
        --condition {params.cond} \
        {params.label} \
        --files {params.files} \
        --replicates {params.replicates} \
        --threshold {params.thresh} \
        --outdir {params.outdir} > {log}
        """


rule statistic_combine_oligo_correlation:
    conda:
        "../../envs/default.yaml"
    input:
        correlation=lambda wc: expand(
            "results/experiments/{{project}}/stats/assigned_counts/{{assignment}}/{{config}}/{condition}_correlation.tsv",
            condition=getConditions(wc.project),
        ),
        correlation_thresh=lambda wc: expand(
            "results/experiments/{{project}}/stats/assigned_counts/{{assignment}}/{{config}}/{condition}_correlation_minThreshold.tsv",
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/experiments/{project}/stats/statistic_oligo_correlation_merged_{assignment}_{config}.tsv",
            caption="../../report/oligo_correlation.rst",
            category="{project}",
            subcategory="Oligos",
        ),
    params:
        thresh=lambda wc: config["experiments"][wc.project]["configs"][wc.config][
            "bc_threshold"
        ],
    log:
        "logs/experiments/{project}/stats/statistic_combine_oligo_correlation.{assignment}_{config}.log",
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
        ) > {output}
        """
