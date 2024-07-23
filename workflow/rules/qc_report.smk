import os

# rule qc_report:
#     input:
#         getOutputProject_helper(
#             [
#                 "results/experiments/{project}/qc_report/qc_report.html",
#             ]
#         ),
#        expand("results/assignment/{assignment}/qc_report.html", assignment=config['assignments'])  

rule qc_report_assoc:
    input: 
        quarto_script = getScript("report/qc_report_assoc.qmd"),
        design_file = lambda wc: config["assignments"][wc.assignment]["reference"], 
        statistic_filter="results/assignment/{assignment}/statistic/assigned_counts.{assignment_config}.tsv",
        statistic_all="results/assignment/{assignment}/statistic/total_counts.tsv",
        plot="results/assignment/{assignment}/statistic/assignment.{assignment_config}.png",
    output: 
        assi_file = "results/assignment/{assignment}/qc_report.{assignment_config}.html",
        quarto_file = temp("results/assignment/{assignment}/qc_report.{assignment_config}.qmd"),
    conda:
        "../envs/quarto.yaml",  
    params:
        bc_length=lambda wc: config["assignments"][wc.assignment]["bc_length"],
        fw = lambda wc: ";".join(list(config["assignments"][wc.assignment]["FW"])),
        rev = lambda wc: ";".join(list(config["assignments"][wc.assignment]["REV"])),
        bc = lambda wc: [";".join(list(config["assignments"][wc.assignment][key])) for key in ["BC", "linker", "linker_length"] if key in config["assignments"][wc.assignment]],
        workdir = os.getcwd(),
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
        quarto_script = getScript("report/qc_report_count.qmd"),
        dna_oligo_coor_min_thre_plot = "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_DNA_pairwise_minThreshold.png",
        rna_oligo_coor_min_thre_plot = "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_RNA_pairwise_minThreshold.png",
        rna_oligo_coor_plot = "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_RNA_pairwise.png",
        dna_oligo_coor_plot = "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_DNA_pairwise.png",
        ratio_oligo_coor_plot = "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_Ratio_pairwise.png",
        ratio_oligo_min_thre_plot = "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_Ratio_pairwise_minThreshold.png",
        statistics_all = "results/experiments/{project}/statistic/statistic_assigned_counts_merged_{assignment}_{config}.tsv",
        per_bar_code_dna = "results/experiments/{project}/statistic/barcode/counts/{condition}_{config}_{type}_perBarcode.png",
        bc_coor_dna = "results/experiments/{project}/statistic/barcode/{raw_or_assigned}/{condition}_{config}_barcode_DNA_pairwise.png",
        bc_coor_rna = "results/experiments/{project}/statistic/barcode/{raw_or_assigned}/{condition}_{config}_barcode_RNA_pairwise.png",
        bc_coor_ratio = "results/experiments/{project}/statistic/barcode/{raw_or_assigned}/{condition}_{config}_barcode_Ratio_pairwise.png",
        # TODO remove bc_coor
        # TODO we need coorelations after assignment. Hide the following files.
        #  rna_oligo_coor_plot = "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_RNA_pairwise.png",
        # dna_oligo_coor_plot = "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_DNA_pairwise.png",
        # ratio_oligo_coor_plot = "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}_Ratio_pairwise.png",
        # TODO also add the minimum threashold value.
        # TODO add some explanation from the documentation about the headers in the table.
        # TODO Total oligos seprately just one 
        # TODO remove extra decimal points.
        # TODO add single table as well.
        # DNA_pearson	RNA_pearson	Ratio_pearson Remove the columns from this table and add this table statistic_oligo_correlation_merged_fromFile_default.
        # TODO add a lightening system. Warning pearson coorelation between replcate 1 and 2 is log.
        # TODO Later, after discussion with Max you can get multiple files for the pngs expanding {condition}.


    output:  
        count_file = "results/experiments/{project}/qc_report/qc_report.{assignment}.{config}.{type}.{condition}.{raw_or_assigned}.html",
        quarto_file = temp("results/experiments/{project}/qc_report/qc_report.{assignment}.{config}.{type}.{condition}.{raw_or_assigned}.qmd"),
    conda:
        "../envs/quarto.yaml",   
    params:
        condition = lambda wildcards: getConditions(wildcards.project),
        workdir = os.getcwd(),

    shell:
        """
        cp config.yml results/experiments/{wildcards.project}/qc_report/config.yml;
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
        -P statistics_all:{input.statistics_all} \
        -P per_bar_code_dna:{input.per_bar_code_dna} \
        -P bc_coor_dna:{input.bc_coor_dna} \
        -P bc_coor_rna:{input.bc_coor_rna} \
        -P bc_coor_ratio:{input.bc_coor_ratio} \
        -P workdir:{params.workdir}
        rm config.yml
        """