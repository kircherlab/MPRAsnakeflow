rule qc_report:
    input:
        getOutputProject_helper(
            [
                "results/experiments/{project}/qc_report/qc_report.html",
            ]
        ),
        expand("results/assignment/{assignment}/qc_report.html", assignment=config['assignments'])  

rule qc_report_assoc:
    input: 
        quarto_script = getScript("report/qc_report_assoc.qmd"),

    output: 
        assi_file = "results/assignment/{assignment}/qc_report.html"
         
    conda:
        "../envs/quarto.yaml",  
    params:
        bc_length=lambda wc: config["assignments"][wc.assignment]["bc_length"],
        fw = lambda wc: config["assignments"][wc.assignment]["FW"],
        rev = lambda wc: config["assignments"][wc.assignment]["REV"],
        bc = lambda wc: config["assignments"][wc.assignment]["BC"],
        reference = lambda wc: config["assignments"][wc.assignment]["reference"],
        configs = lambda wc: list(config["assignments"][wc.assignment]["configs"].keys())[0],
        
    shell:
        """
        cd results/assignment/{wildcards.assignment}/
        cp {input.quarto_script} qc_report_assoc.qmd
        quarto render qc_report_assoc.qmd --output qc_report.html \
        -P assignment:{wildcards.assignment} \
        -P bc_length:{params.bc_length} \
        -P fw:{params.fw} \
        -P rev:{params.rev} \
        -P bc:{params.bc} \
        -P reference:{params.reference} \
        -P configs:{params.configs}
        rm qc_report_assoc.qmd
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