####################################
# Statistic of barcodes and oligos #
####################################


##################
## subworkflows ##
##################


# statistic on BC counts (not assigned)
include: "statistic/counts.smk"
# statistic on assigned counts (BCs and oligos)
include: "statistic/assigned_counts.smk"
# statistic on correlation of BCs/coligos
include: "statistic/correlation.smk"
# BC overlap between replicates
include: "statistic/bc_overlap.smk"


rule experiment_statistic_quality_metric:
    """
Quality metrics of the assignment run
"""
    input:
        barcode="results/experiments/{project}/reporter_experiment.barcode.{condition}.{assignment}.{config}.all.tsv.gz",
        assignment="results/experiments/{project}/assignment/{assignment}.tsv.gz",
        script=getScript("quality_metrics.py"),
    output:
        "results/experiments/{project}/qc_metrics.{condition}.{assignment}.{config}.json",
    log:
        "results/logs/experiment/statistic_quality_metric.{project}.{condition}.{assignment}.{config}.log",
    conda:
        getCondaEnv("mpralib.yaml")
    params:
        bc_threshold=lambda wc: config["experiments"][wc.project]["configs"][wc.config]["filter"]["bc_threshold"],
    shell:
        """
        python {input.script} experiment \
        --assignment {input.assignment} --barcode {input.barcode} \
        --bc-threshold {params.bc_threshold} \
        --output {output} &> {log}
        """
