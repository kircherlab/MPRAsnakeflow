rule render_quarto_from_smk:
    output:
        html_file="results/qc_report/output_from_smk_file.html"
    shell:
        "touch results/qc_report/output_from_smk_file.html"