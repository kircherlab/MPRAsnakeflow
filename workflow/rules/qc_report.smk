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
        condition = getOutputProjectConditionConfigType_helper(["{condition}"]),
        config = getOutputProjectConditionConfigType_helper(["{config}"]),
        type = getOutputProjectConditionConfigType_helper(["{type}"]),
        
    shell:
        """ 
        echo "Available wildcards: {wildcards}"
        echo "testing condition: "
        echo {params.condition}
        echo "testing config: "
        echo {params.config}
        echo "testing type: "
        echo {params.type}
        cd results/experiments/{wildcards.project}/qc_report
        cp {input.quarto_script} qc_report.qmd
        quarto render qc_report.qmd --output qc_report.html \
        -P image_url:{params.perbarcode_dna} \
        -P condition:{params.condition} \
        -P config:{params.config}
        rm qc_report.qmd
        """