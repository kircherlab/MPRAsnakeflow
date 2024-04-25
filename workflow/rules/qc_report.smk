
## Function definations.
def getProjects():
    if "experiments" in config:
        return list(config["experiments"].keys())
    else:
        return []

def getOutputProject_helper(files, betweenReplicates=False):
    """
    Inserts {project}, {config} from config into given file.
    When betweenReplicates is True skips projects without replicates in one condition.
    """
    output = []
    projects = getProjects()
    for project in projects:
        if not betweenReplicates or hasReplicates(project):
            for file in files:
                output += expand(
                    file,
                    project=project,
                )
    return output


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
    shell:
        """ 
        cd results/experiments/{wildcards.project}/qc_report
        cp {input.quarto_script} qc_report.qmd
        quarto render qc_report.qmd --output qc_report.html -P image_url:HEPG2.png
        rm qc_report.qmd
        """