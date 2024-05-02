rule qc_report:
    input:
        getOutputProject_helper(
            [
                "results/experiments/{project}/qc_report/qc_report.html",
            ]
        ),


rule report_generator:
    input: 
        quarto_script = getScript("report/qc_report.qmd"),
    output:
        "results/experiments/{project}/qc_report/qc_report.html",
    conda:
        "../envs/quarto.yaml",   
    params:
        condition = lambda wildcards: getConditions(wildcards.project),
        # config = getOutputProjectConditionConfigType_helper(["{config}"]),
        # types = getOutputProjectConditionConfigType_helper(["{type}"]),
        
    shell:
        """ 
        cp config.yml results/experiments/{wildcards.project}/qc_report/config.yml
        cd results/experiments/{wildcards.project}/qc_report
        cp {input.quarto_script} qc_report.qmd
        quarto render qc_report.qmd --output qc_report.html \
        -P condition:{params.condition} \
        -P project:{wildcards.project} \
        --execute-params config.yml 
        rm qc_report.qmd
        """


        #  -P image_url:{params.perbarcode_dna} \