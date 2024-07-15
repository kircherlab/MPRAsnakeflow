import os


rule qc_report:
    input:
        getOutputProject_helper(
            [
                "results/experiments/{project}/qc_report/qc_report.html",
            ]
        ),
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
        fw = lambda wc: config["assignments"][wc.assignment]["FW"],
        rev = lambda wc: config["assignments"][wc.assignment]["REV"],
        bc = lambda wc: config["assignments"][wc.assignment]["BC"],
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
        quarto_script = getScript("report/qc_report.qmd"),
    output:
        "results/experiments/{project}/qc_report/qc_report.html",
    conda:
        "../envs/quarto.yaml",   
    params:
        condition = lambda wildcards: getConditions(wildcards.project),

    shell:
        """
        cp config.yml results/experiments/{wildcards.project}/qc_report/config.yml
        cd results/experiments/{wildcards.project}/qc_report
        cp {input.quarto_script} qc_report.qmd
        quarto render qc_report.qmd --output qc_report.html \
        -P condition:{params.condition} \
        -P project:{wildcards.project}
        rm qc_report.qmd
        rm config.yml
        """
