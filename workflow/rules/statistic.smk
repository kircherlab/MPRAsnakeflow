# Statistic of barcodes and oligos

#################################
## Count statistic of barcodes ##
#################################


# count Reads, Barcodes per UMI, Barcodes and Unique UMIs
rule statistic_counts:
    input:
        "results/{project}/counts/{condition}_{replicate}_{type}_{countType}_counts.tsv.gz",
    output:
        "results/{project}/stats/counts/{condition}_{replicate}_{type}_{countType}_counts.tsv.gz",
    params:
        cond="{condition}",
        rep="{replicate}",
        type="{type}",
    shell:
        """
        paste <( echo "{params.cond}") <( echo "{params.rep}") <( echo "{params.type}") \
        <( 
            zcat {input} | \
            awk -v OFS='\\t' 'BEGIN{{
                pbar="NA"
            }}{{ 
                count += $NF; umi_sum+=$3; if (pbar != $1) {{ barcodes+=1 }}; pbar=$1 
            }}END{{ 
                if (NR > 0) {{
                    print umi_sum/NR,count,NR,barcodes
                }} else {{
                    print 0,0,0,0
                }}
            }}' 
        ) \
        <( 
            zcat {input} | cut -f 2 | sort -u | wc -l 
        ) | \
        gzip -c > {output}
        """


# get all counts of experiment (rule statistic_counts)
def getCountStats(wc):
    exp = getExperiments(wc.project)
    output = []
    for index, row in exp.iterrows():
        output += expand(
            "results/{project}/stats/counts/{condition}_{replicate}_{type}_{countType}_counts.tsv.gz",
            project=wc.project,
            condition=row["Condition"],
            replicate=row["Replicate"],
            type=["DNA", "RNA"],
            countType=wc.countType,
        )
    return output


# concat DNA, RNA-counts (rule statistic_counts) for all experiments, and replicates
rule count_stats_merge:
    input:
        getCountStats,
    output:
        "results/{project}/stats/counts/count_{countType}.tsv",
    shell:
        """
        zcat {input} | sort -k1,1 -k3,3 -k2,2 > {output}
        """


# Statistic of barcodes shared between RNA&DNA per condition and replicate
rule statistic_BC_in_RNA_DNA:
    input:
        dna="results/{project}/counts/{condition}_{replicate}_DNA_{countType}_counts.tsv.gz",
        rna="results/{project}/counts/{condition}_{replicate}_RNA_{countType}_counts.tsv.gz",
    output:
        "results/{project}/stats/counts/{condition}_{replicate}_{countType}_BC_in_RNA_DNA.tsv.gz",
    params:
        cond="{condition}",
        rep="{replicate}",
    shell:
        """
        paste <( echo "{params.cond}") <( echo "{params.rep}") \
        <( join <( zcat {input.dna} | cut -f 1 | sort | uniq ) \
        <( zcat {input.rna} | cut -f 1 | sort | uniq ) | wc -l ) | \
        gzip -c > {output}
        """


# get all barcodes of experiment (rule statistic_BC_in_RNA_DNA)
def getBCinRNADNAStats(wc):
    exp = getExperiments(wc.project)
    output = []
    for index, row in exp.iterrows():
        output += expand(
            "results/{project}/stats/counts/{condition}_{replicate}_{countType}_BC_in_RNA_DNA.tsv.gz",
            project=wc.project,
            condition=row["Condition"],
            replicate=row["Replicate"],
            countType=wc.countType,
        )
    return output


# concat shared barcodes (rule statistic_BC_in_RNA_DNA) for all experiments, and replicates
rule stats_BC_in_RNA_DNA_merge:
    input:
        getBCinRNADNAStats,
    output:
        "results/{project}/stats/counts/BC_in_RNA_DNA_{countType}.tsv",
    shell:
        """
        zcat {input} | sort -k1,1 -k2,2 > {output}
        """


# making final count statistics
rule stats_final:
    conda:
        "../envs/mpraflow_r.yaml"
    input:
        counts="results/{project}/stats/counts/count_{countType}.tsv",
        shared="results/{project}/stats/counts/BC_in_RNA_DNA_{countType}.tsv",
    output:
        "results/{project}/stats/statistic_count_{countType}.tsv",
    shell:
        """
        Rscript workflow/scripts/count/combine_count_stats.R --count {input.counts} --shared {input.shared} --output {output}
        """


#################################
## count 10 most frequent UMIs ##
#################################


# count frequent UMIs per condition, replicate and DNA/RNA
rule frequent_umis:
    input:
        "results/{project}/counts/{condition}_{replicate}_{type}_filtered_counts.tsv.gz",
    output:
        freqUMIs=(
            "results/{project}/stats/counts/freqUMIs_{condition}_{replicate}_{type}.txt"
        ),
    shell:
        """
        set +o pipefail;
        zcat {input} | cut -f 2 | sort | uniq -c | sort -nr | head > {output.freqUMIs}
        """


##############################
## Barcode base composition ##
##############################


rule barcode_base_composition:
    conda:
        "../envs/mpraflow_pandas.yaml"
    input:
        counts="results/{project}/counts/{condition}_{replicate}_{type}_final_counts.tsv.gz",
    output:
        bc=temp(
            "results/{project}/counts/{condition}_{replicate}_{type}_final.BC.tsv.gz"
        ),
        stats="results/{project}/stats/counts/BCNucleotideComposition/{condition}_{replicate}_{type}.tsv.gz",
    params:
        name="{condition}_{replicate}_{type}",
    shell:
        """
        zcat {input.counts} | awk '{{print $2}}' | gzip -c > {output.bc};
        python workflow/scripts/count/nucleotideCountPerPosition.py \
        --column 1 \
        --chunksize 100000 \
        --input {output.bc} \
        --output {output.stats}
        """


#############################
## Correlation of Barcodes ##
#############################


# overlap barcodes and counts dna and RNA seperately

# first assign RNA and DNA barcodes seperately to make the statistic for assigned


rule assignBarcodes:
    conda:
        "../envs/mpraflow_py36.yaml"
    input:
        counts="results/{project}/counts/{condition}_{replicate}_{type}_final_counts.tsv.gz",
        association=lambda wc: config[wc.project]["assignments"][wc.assignment],
    output:
        counts="results/{project}/assigned_counts/{assignment}/{condition}_{replicate}_{type}_final_counts.tsv.gz",
        stats="results/{project}/stats/assigned_counts/{assignment}/{condition}_{replicate}_{type}.statistic.tsv.gz",
    params:
        name="{condition}_{replicate}_{type}",
    shell:
        """
        python workflow/scripts/count/merge_BC_and_assignment.py --counts {input.counts} \
        --assignment {input.association} \
        --output {output.counts} \
        --statistic {output.stats} \
        --name {params.name}
        """


rule combine_BC_assignment_stats_helper:
    input:
        stats=lambda wc: expand(
            "results/{{project}}/stats/assigned_counts/{{assignment}}/{{condition}}_{replicate}_{type}.statistic.tsv.gz",
            type=["DNA", "RNA"],
            replicate=getReplicatesOfCondition(wc.project, wc.condition),
        ),
    output:
        temp(
            "results/{project}/stats/assigned_counts/{assignment}/helper.{condition}.statistic.tsv.gz"
        ),
    shell:
        """
        set +o pipefail;
        (
            zcat {input.stats[0]} | head -n 1;
            for i in {input.stats}; do
                zcat $i | tail -n +2
            done;
        ) | gzip -c > {output}
        """


rule combine_BC_assignment_stats:
    input:
        stats=lambda wc: expand(
            "results/{{project}}/stats/assigned_counts/{{assignment}}/helper.{condition}.statistic.tsv.gz",
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/{project}/stats/statistic_assigned_counts_single_{assignment}.tsv",
            caption="../report/assigned_counts_beforeMerge.rst",
            category="{project}",
            subcategory="Assignment",
        ),
    shell:
        """
        set +o pipefail;
        (
            zcat {input.stats[0]} | head -n 1;
            for i in {input.stats}; do
                zcat $i | tail -n +2
            done;
        ) > {output}
        """


###


rule overlapBCs:
    conda:
        "../envs/mpraflow_r.yaml"
    input:
        lambda wc: expand(
            "results/{{project}}/{{raw_or_assigned}}/{{condition}}_{replicate}_{{type}}_final_counts.tsv.gz",
            replicate=getReplicatesOfCondition(wc.project, wc.condition),
        ),
    output:
        "results/{project}/stats/{raw_or_assigned}/overlapBCandCounts_{condition}_{type}.tsv",
    params:
        input=lambda wc: ",".join(
            expand(
                "results/{project}/{raw_or_assigned}/{condition}_{replicate}_{type}_final_counts.tsv.gz",
                project=wc.project,
                condition=wc.condition,
                raw_or_assigned=wc.raw_or_assigned,
                type=wc.type,
                replicate=getReplicatesOfCondition(wc.project, wc.condition),
            )
        ),
        cond="{condition}_{type}",
        replicates=lambda wc: ",".join(
            getReplicatesOfCondition(wc.project, wc.condition)
        ),
    shell:
        """
        Rscript workflow/scripts/count/BCCounts_betweenReplicates.R \
        --outfile {output} \
        --condition {params.cond} \
        --files {params.input} --replicates {params.replicates}
        """


rule combine_overlapBCs_stats_raw:
    input:
        stats=lambda wc: expand(
            "results/{{project}}/stats/counts/overlapBCandCounts_{condition}_{type}.tsv",
            type=["DNA", "RNA"],
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/{project}/stats/statistic_overlapBCs_counts.tsv",
            caption="../report/bc_overlap.rst",
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


rule combine_overlapBCs_stats_assigned:
    input:
        stats=lambda wc: expand(
            "results/{{project}}/stats/assigned_counts/{{assignment}}/overlapBCandCounts_{condition}_{type}.tsv",
            type=["DNA", "RNA"],
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/{project}/stats/statistic_overlapBCs_assigned_counts_{assignment}.tsv",
            caption="../report/bc_overlap_assignment.rst",
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


# get all barcodes of experiment (rule dna_rna_merge_counts_withoutZeros or rule dna_rna_merge_counts_withZeros)
def getMergedCounts(wc):
    exp = getExperiments(wc.project)
    exp = exp[exp.Condition == wc.condition]
    files = []
    replicates = []
    for index, row in exp.iterrows():
        files += expand(
            "results/{project}/{raw_or_assigned}/merged/{mergeType}/{condition}_{replicate}_merged_counts.tsv.gz",
            raw_or_assigned=wc.raw_or_assigned,
            project=wc.project,
            condition=wc.condition,
            replicate=row["Replicate"],
            mergeType=wc.mergeType,
        )
        replicates += str(row["Replicate"])
    return [files, replicates]


rule correlate_BC_counts:
    conda:
        "../envs/mpraflow_r.yaml"
    input:
        lambda wc: getMergedCounts(wc)[0],
    output:
        "results/{project}/stats/barcode/{raw_or_assigned}/{mergeType}/{condition}_barcode_DNA_pairwise.png",
        "results/{project}/stats/barcode/{raw_or_assigned}/{mergeType}/{condition}_barcode_RNA_pairwise.png",
        "results/{project}/stats/barcode/{raw_or_assigned}/{mergeType}/{condition}_barcode_Ratio_pairwise.png",
        "results/{project}/stats/barcode/{raw_or_assigned}/{mergeType}/{condition}_barcode_correlation.tsv",
        "results/{project}/stats/barcode/{raw_or_assigned}/{mergeType}/{condition}_DNA_perBarcode.png",
        "results/{project}/stats/barcode/{raw_or_assigned}/{mergeType}/{condition}_RNA_perBarcode.png",
    params:
        replicates=lambda wc: ",".join(getMergedCounts(wc)[1]),
        cond="{condition}",
        outdir="results/{project}/stats/barcode/{raw_or_assigned}/{mergeType}",
        input=lambda wc: ",".join(getMergedCounts(wc)[0]),
    shell:
        """
        Rscript workflow/scripts/count/plot_perBCCounts_correlation.R \
        --outdir {params.outdir} \
        --condition {params.cond} \
        --files {params.input} --replicates {params.replicates}
        """


rule combine_bc_correlation_raw:
    conda:
        "../envs/mpraflow_py36.yaml"
    input:
        lambda wc: expand(
            "results/{{project}}/stats/barcode/counts/{{mergeType}}/{condition}_barcode_correlation.tsv",
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/{project}/stats/statistic_bc_correlation_merged_{mergeType}.tsv",
            caption="../report/bc_correlation.rst",
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


rule combine_bc_correlation_assigned:
    conda:
        "../envs/mpraflow_py36.yaml"
    input:
        lambda wc: expand(
            "results/{{project}}/stats/barcode/assigned_counts/{{assignment}}/{{mergeType}}/{condition}_barcode_correlation.tsv",
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/{project}/stats/statistic_assigned_bc_correlation_merged_{mergeType}_{assignment}.tsv",
            caption="../report/bc_correlation_assigned.rst",
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


##########


def getAssignedCountsStatistic(wc):
    exp = getExperiments(wc.project)
    exp = exp[exp.Condition == wc.condition]
    output = []
    for index, row in exp.iterrows():
        output += [
            "--statistic %s results/%s/stats/assigned_counts/%s/%s/%s_%s_merged_assigned_counts.statistic.tsv.gz"
            % (
                str(row["Replicate"]),
                wc.project,
                wc.assignment,
                wc.config,
                wc.condition,
                str(row["Replicate"]),
            )
        ]
    return output


rule combine_stats_dna_rna_merge:
    conda:
        "../envs/mpraflow_py36.yaml"
    input:
        lambda wc: expand(
            "results/{{project}}/stats/assigned_counts/{{assignment}}/{{config}}/{{condition}}_{replicate}_merged_assigned_counts.statistic.tsv.gz",
            replicate=getReplicatesOfCondition(wc.project, wc.condition),
        ),
    output:
        "results/{project}/stats/assigned_counts/{assignment}/{config}/combined/{condition}_merged_assigned_counts.statistic.tsv.gz",
    params:
        cond="{condition}",
        statistic=lambda wc: " ".join(getAssignedCountsStatistic(wc)),
    shell:
        """
        python workflow/scripts/count/merge_statistic_tables.py \
        --condition {params.cond} \
        {params.statistic} \
        --output {output}
        """


rule combine_stats_dna_rna_merge_all:
    conda:
        "../envs/mpraflow_py36.yaml"
    input:
        lambda wc: expand(
            "results/{{project}}/stats/assigned_counts/{{assignment}}/{{config}}/combined/{condition}_merged_assigned_counts.statistic.tsv.gz",
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/{project}/stats/statistic_assigned_counts_merged_{assignment}_{config}.tsv",
            caption="../report/assigned_counts_afterMerge.rst",
            category="{project}",
            subcategory="Assignment",
        ),
    shell:
        """
        set +o pipefail;
        (
            zcat {input[0]} | head -n 1;
            for i in {input}; do
                zcat $i | tail -n +2
            done
        ) > {output}
        """


##############


rule calc_correlations:
    conda:
        "../envs/mpraflow_r.yaml"
    input:
        counts=lambda wc: expand(
            "results/{{project}}/assigned_counts/{{assignment}}/{{config}}/{{condition}}_{replicate}_merged_assigned_counts.tsv.gz",
            replicate=getReplicatesOfCondition(wc.project, wc.condition),
        ),
        label=(
            lambda wc: config[wc.project]["label_file"]
            if "label_file" in config[wc.project]
            else []
        ),
    output:
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_all_barcodesPerInsert_box.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_all_barcodesPerInsert_box_minThreshold.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_DNA_pairwise.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_DNA_pairwise_minThreshold.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_group_barcodesPerInsert_box.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_group_barcodesPerInsert_box_minThreshold.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_Ratio_pairwise.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_Ratio_pairwise_minThreshold.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_RNA_pairwise.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_RNA_pairwise_minThreshold.png",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_correlation.tsv",
        "results/{project}/stats/assigned_counts/{assignment}/{config}/{condition}_correlation_minThreshold.tsv",
    params:
        cond="{condition}",
        files=lambda wc: ",".join(
            expand(
                "results/{project}/assigned_counts/{assignment}/{config}/{condition}_{replicate}_merged_assigned_counts.tsv.gz",
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
        thresh=lambda wc: config[wc.project]["configs"][wc.config]["bc_threshhold"],
        outdir="results/{project}/stats/assigned_counts/{assignment}/{config}",
        label=(
            lambda wc: "--label %s" % config[wc.project]["label_file"]
            if "label_file" in config[wc.project]
            else ""
        ),
    shell:
        """
        Rscript workflow/scripts/count/plot_perInsertCounts_correlation.R \
        --condition {params.cond} \
        {params.label} \
        --files {params.files} \
        --replicates {params.replicates} \
        --threshold {params.thresh} \
        --outdir {params.outdir}
        """


rule combine_oligo_correlation:
    conda:
        "../envs/mpraflow_py36.yaml"
    input:
        correlation=lambda wc: expand(
            "results/{{project}}/stats/assigned_counts/{{assignment}}/{{config}}/{condition}_correlation.tsv",
            condition=getConditions(wc.project),
        ),
        correlation_thresh=lambda wc: expand(
            "results/{{project}}/stats/assigned_counts/{{assignment}}/{{config}}/{condition}_correlation_minThreshold.tsv",
            condition=getConditions(wc.project),
        ),
    output:
        report(
            "results/{project}/stats/statistic_oligo_correlation_merged_{assignment}_{config}.tsv",
            caption="../report/oligo_correlation.rst",
            category="{project}",
            subcategory="Oligos",
        ),
    params:
        thresh=lambda wc: config[wc.project]["configs"][wc.config]["bc_threshhold"],
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
