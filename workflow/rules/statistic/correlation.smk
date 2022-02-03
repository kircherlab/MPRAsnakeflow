################################
### Correlation of BC Counts ###
################################

# get all barcodes of experiment (rule dna_rna_merge_counts_withoutZeros or rule dna_rna_merge_counts_withZeros)
def getMergedCounts(project, raw_or_assigned, condition, mergeType, sampling):
    exp = getExperiments(project)
    exp = exp[exp.Condition == condition]
    files = []
    replicates = []
    for index, row in exp.iterrows():
        files += expand(
            "results/{project}/{raw_or_assigned}/merged/{mergeType}/{condition}_{replicate}_merged_counts_{sampling}.tsv.gz",
            raw_or_assigned=raw_or_assigned,
            project=project,
            condition=condition,
            replicate=row["Replicate"],
            sampling=sampling,
            mergeType=mergeType,
        )
        replicates += str(row["Replicate"])
    return [files, replicates]


rule statistic_correlate_BC_counts:
    conda:
        "../../envs/mpraflow_r.yaml"
    input:
        lambda wc: getMergedCounts(
            wc.project, wc.raw_or_assigned, wc.condition, wc.mergeType, wc.sampling
        )[0],
    output:
        "results/{project}/stats/barcode/{raw_or_assigned}/{mergeType}/{condition}_{sampling}_barcode_DNA_pairwise.png",
        "results/{project}/stats/barcode/{raw_or_assigned}/{mergeType}/{condition}_{sampling}_barcode_RNA_pairwise.png",
        "results/{project}/stats/barcode/{raw_or_assigned}/{mergeType}/{condition}_{sampling}_barcode_Ratio_pairwise.png",
        "results/{project}/stats/barcode/{raw_or_assigned}/{mergeType}/{condition}_{sampling}_barcode_correlation.tsv",
        "results/{project}/stats/barcode/{raw_or_assigned}/{mergeType}/{condition}_{sampling}_DNA_perBarcode.png",
        "results/{project}/stats/barcode/{raw_or_assigned}/{mergeType}/{condition}_{sampling}_RNA_perBarcode.png",
    params:
        replicates=lambda wc: ",".join(
            getMergedCounts(
                wc.project, wc.raw_or_assigned, wc.condition, wc.mergeType, wc.sampling
            )[1]
        ),
        cond="{condition}",
        outdir="results/{project}/stats/barcode/{raw_or_assigned}/{mergeType}",
        input=lambda wc: ",".join(
            getMergedCounts(
                wc.project, wc.raw_or_assigned, wc.condition, wc.mergeType
            )[0]
        ),
    shell:
        """
        Rscript {SCRIPTS_DIR}/count/plot_perBCCounts_correlation.R \
        --outdir {params.outdir} \
        --condition {params.cond} \
        --files {params.input} --replicates {params.replicates}
        """


rule statistic_combine_bc_correlation_raw:
    conda:
        "../../envs/mpraflow_py36.yaml"
    input:
        lambda wc: expand(
            "results/{{project}}/stats/barcode/counts/{{mergeType}}/{condition}_{{sampling}}_barcode_correlation.tsv",
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/{project}/stats/statistic_bc_correlation_merged_{mergeType}_{sampling}.tsv",
            caption="../../report/bc_correlation.rst",
            category="{project}",
            subcategory="Barcodes",
        ),
    shell:
        """
        (
        cat {input[0]} | head -n 1;
        for i in {input}; do
            cat $i | tail -n +2;
        done;
        ) > {output}
        """


rule statistic_combine_bc_correlation_assigned:
    conda:
        "../../envs/mpraflow_py36.yaml"
    input:
        lambda wc: expand(
            "results/{{project}}/stats/barcode/assigned_counts/{{assignment}}/{{mergeType}}/{condition}_{{sampling}}_barcode_correlation.tsv",
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/{project}/stats/statistic_assigned_bc_correlation_merged_{mergeType}_{assignment}_{sampling}.tsv",
            caption="../../report/bc_correlation_assigned.rst",
            category="{project}",
            subcategory="Barcodes",
        ),
    shell:
        """
        (
        cat {input[0]} | head -n 1;
        for i in {input}; do
            cat $i | tail -n +2;
        done;
        ) > {output}
        """


#############################
### Correlation of Oligos ###
#############################


rule statistic_calc_correlations:
    conda:
        "../../envs/mpraflow_r.yaml"
    input:
        counts=lambda wc: expand(
            "results/{{project}}/assigned_counts/{{assignment}}/{{config}}/{{condition}}_{replicate}_merged_assigned_counts_{{sampling}}.tsv.gz",
            replicate=getReplicatesOfCondition(wc.project, wc.condition),
        ),
        label=(
            lambda wc: config[wc.project]["label_file"]
            if "label_file" in config[wc.project]
            else []
        ),
    output:
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_{sampling}_all_barcodesPerInsert_box.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_{sampling}_all_barcodesPerInsert_box_minThreshold.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_{sampling}_DNA_pairwise.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_{sampling}_DNA_pairwise_minThreshold.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_{sampling}_group_barcodesPerInsert_box.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_{sampling}_group_barcodesPerInsert_box_minThreshold.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_{sampling}_Ratio_pairwise.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_{sampling}_Ratio_pairwise_minThreshold.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_{sampling}_RNA_pairwise.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_{sampling}_RNA_pairwise_minThreshold.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_{sampling}_correlation.tsv",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_{sampling}_correlation_minThreshold.tsv",
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
        outdir="results/{project}/stats/assigned_counts/{assignment}/{config}",
        label=(
            lambda wc: "--label %s" % config[wc.project]["label_file"]
            if "label_file" in config[wc.project]
            else ""
        ),
    shell:
        """
        Rscript {SCRIPTS_DIR}/count/plot_perInsertCounts_correlation.R \
        --condition {params.cond} \
        {params.label} \
        --files {params.files} \
        --replicates {params.replicates} \
        --threshold {params.thresh} \
        --outdir {params.outdir}
        """


rule statistic_combine_oligo_correlation:
    conda:
        "../../envs/mpraflow_py36.yaml"
    input:
        correlation=lambda wc: expand(
            "results/{{project}}/stats/assigned_counts/{{assignment}}/{{config}}/{condition}_{{sampling}}_correlation.tsv",
            condition=getConditions(wc.project),
        ),
        correlation_thresh=lambda wc: expand(
            "results/{{project}}/stats/assigned_counts/{{assignment}}/{{config}}/{condition}_{{sampling}}_correlation_minThreshold.tsv",
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/{project}/stats/statistic_oligo_correlation_merged_{assignment}_{config}_{sampling}.tsv",
            caption="../../report/oligo_correlation.rst",
            category="{project}",
            subcategory="Oligos",
        ),
    params:
        thresh=lambda wc: config[wc.project]["configs"][wc.config]["bc_threshold"],
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
