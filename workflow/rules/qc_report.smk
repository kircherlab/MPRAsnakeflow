import os


rule qc_report_assoc:
    input:
        quarto_script=getScript("report/qc_report_assoc.qmd"),
        design_file=lambda wc: config["assignments"][wc.assignment]["design_file"],
        statistic_filter="results/assignment/{assignment}/statistic/assigned_counts.{assignment_config}.tsv",
        statistic_all="results/assignment/{assignment}/statistic/total_counts.tsv",
        plot="results/assignment/{assignment}/statistic/assignment.{assignment_config}.png",
    output:
        assi_file="results/assignment/{assignment}/qc_report.{assignment_config}.html",
        quarto_file=temp(
            "results/assignment/{assignment}/qc_report.{assignment_config}.qmd"
        ),
    conda:
        "../envs/quarto.yaml"
    params:
        bc_length=lambda wc: config["assignments"][wc.assignment]["bc_length"],
        fw=lambda wc: ";".join(list(config["assignments"][wc.assignment]["FW"])),
        rev=lambda wc: ";".join(list(config["assignments"][wc.assignment]["REV"])),
        bc=lambda wc: [
            ";".join(list(config["assignments"][wc.assignment][key]))
            for key in ["BC", "linker", "linker_length"]
            if key in config["assignments"][wc.assignment]
        ],
        workdir=os.getcwd(),
    shell:
        """
        cp {input.quarto_script} {output.quarto_file};
        cd `dirname {output.quarto_file}`;
        quarto render `basename {output.quarto_file}` --output `basename {output.assi_file}` \
        -P assignment:{wildcards.assignment} \
        -P bc_length:{params.bc_length} \
        -P fw:{params.fw} \
        -P rev:{params.rev} \
        -P bc:{params.bc} \
        -P workdir:{params.workdir} \
        -P design_file:{input.design_file} \
        -P configs:{wildcards.assignment_config} \
        -P plot_file:{input.plot} \
        -P statistic_filter_file:{input.statistic_filter} \
        -P statistic_all_file:{input.statistic_all}
        """


rule qc_report_count:
    input:
        quarto_script=getScript("report/qc_report_count.qmd"),
        dna_oligo_coor_min_thre_plot="results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_DNA_pairwise_minThreshold.png",
        rna_oligo_coor_min_thre_plot="results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_RNA_pairwise_minThreshold.png",
        rna_oligo_coor_plot="results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_RNA_pairwise.png",
        dna_oligo_coor_plot="results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_DNA_pairwise.png",
        ratio_oligo_coor_plot="results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_Ratio_pairwise.png",
        ratio_oligo_min_thre_plot="results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_Ratio_pairwise_minThreshold.png",
        statistics_all_merged="results/experiments/{project}/statistic/statistic_assigned_counts_merged_{assignment}_{config}.tsv",
        statistics_all_single="results/experiments/{project}/statistic/statistic_assigned_counts_single_{assignment}_{config}.tsv",
        statistics_all_oligo_cor_all="results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_correlation.tsv",
        statistics_all_oligo_cor_thresh="results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_correlation_minThreshold.tsv",
        counts_per_oligo_dna="results/experiments/{project}/statistic/barcode/assigned_counts/{assignment}/{condition}_{config}_DNA_perBarcode.png",
        counts_per_oligo_rna="results/experiments/{project}/statistic/barcode/assigned_counts/{assignment}/{condition}_{config}_RNA_perBarcode.png",
        activity_thresh="results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_group_barcodesPerInsert_box.png",
        activity_all="results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_group_barcodesPerInsert_box_minThreshold.png",
        # TODO Later, after discussion with Max you can get multiple files for the pngs expanding {condition}.
    output:
        count_file="results/experiments/{project}/qc_report.{condition}.{assignment}.{config}.html",
        quarto_file=temp(
            "results/experiments/{project}/qc_report.{condition}.{assignment}.{config}.qmd"
        ),
    conda:
        "../envs/quarto.yaml"
    params:
        condition=lambda wildcards: getConditions(wildcards.project),
        workdir=os.getcwd(),
        thresh=lambda wildcards: str(
            config["experiments"][wildcards.project]["configs"][wildcards.config][
                "filter"
            ]["bc_threshold"]
        ),
    shell:
        """
        cp {input.quarto_script} {output.quarto_file};
        cd `dirname {output.quarto_file}`;
        quarto render `basename {output.quarto_file}` --output `basename {output.count_file}` \
        -P assignment:{wildcards.assignment} \
        -P project:{wildcards.project} \
        -P dna_oligo_coor_min_thre_plot:{input.dna_oligo_coor_min_thre_plot} \
        -P rna_oligo_coor_min_thre_plot:{input.rna_oligo_coor_min_thre_plot} \
        -P dna_oligo_coor_plot:{input.dna_oligo_coor_plot} \
        -P rna_oligo_coor_plot:{input.rna_oligo_coor_plot} \
        -P ratio_oligo_coor_plot:{input.ratio_oligo_coor_plot} \
        -P ratio_oligo_min_thre_plot:{input.ratio_oligo_min_thre_plot} \
        -P statistics_all_merged:{input.statistics_all_merged} \
        -P counts_per_oligo_dna:{input.counts_per_oligo_dna} \
        -P counts_per_oligo_rna:{input.counts_per_oligo_rna} \
        -P statistics_all_single:{input.statistics_all_single} \
        -P activity_all:{input.activity_all} \
        -P activity_thresh:{input.activity_thresh} \
        -P statistics_all_oligo_cor_all:{input.statistics_all_oligo_cor_all} \
        -P statistics_all_oligo_cor_thresh:{input.statistics_all_oligo_cor_thresh} \
        -P thresh:{params.thresh} \
        -P workdir:{params.workdir}
        """