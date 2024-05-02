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
        my_param = "some_value",  
        perbarcode_dna = [s for s in
            
            getOutputProjectConditionConfigType_helper(
            [
                "/statistic/barcode/counts/{condition}_{config}_{type}_perBarcode.png"
                ]
            )
            if "DNA" in s],
        perbarcode = getOutputProjectConditionConfigType_helper(["/statistic/barcode/counts/{condition}_{config}_{type}_perBarcode.png"]),
        # condition = getOutputProjectConditionConfigType_helper(["{condition}"]),
        # condition = getConditions("{wildcards.project}"),
        condition = lambda wildcards: getConditions(wildcards.project),

        config = getOutputProjectConditionConfigType_helper(["{config}"]),
        types = getOutputProjectConditionConfigType_helper(["{type}"]),
        
    shell:
        """ 
        echo "Available wildcards: {wildcards}"
        echo "testing condition: {params.condition}"
        echo "testing config: {params.config}"
        echo "checking types: {params.types}"
        ls
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