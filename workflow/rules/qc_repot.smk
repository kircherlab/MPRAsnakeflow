rule render_quarto:
    input:
        # qmd_file="minimal.qmd"
        qmd_file=getScript("report/minimal.qmd"),

        # workflow/scripts/report/minimal.qmd
    output:
        html_file="results/experiments/{project}/output.html"
        
    conda:
        "../envs/quarto.yaml"
    shell:
        "quarto render {input.qmd_file} --output {output.html_file}"

