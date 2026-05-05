##################################
### Assign BCs and afterwards ###
##################################


rule experiment_assigned_counts_filterAssignment:
    """
Use only unique assignments and do sampling if needed.
"""
    input:
        assignment=lambda wc: getAssignmentFile(wc.project, wc.assignment),
        script=getScript("count/samplerer_assignment.py"),
    output:
        "results/experiments/{project}/assignment/{assignment}.tsv.gz",
    log:
        temp("results/logs/experiment/assigned_counts/filterAssignment.{project}.{assignment}.log"),
    conda:
        getCondaEnv("python3.yaml")
    params:
        samplingprop=lambda wc: assignedCounts_getAssignmentSamplingConfig(wc.project, wc.assignment, "prop"),
        samplingtotal=lambda wc: assignedCounts_getAssignmentSamplingConfig(wc.project, wc.assignment, "total"),
        seed=lambda wc: assignedCounts_getAssignmentSamplingConfig(wc.project, wc.assignment, "seed"),
    shell:
        """
        python {input.script} \
        --input {input.assignment} \
        {params.samplingprop} \
        {params.samplingtotal} \
        {params.seed} \
        --output {output} &> {log}
        """


rule experiment_assigned_counts_createAssignmentPickleFile:
    input:
        files="results/experiments/{project}/assignment/{assignment}.tsv.gz",
        script=getScript("count/create_pickle.py"),
    output:
        "results/experiments/{project}/assignment/{assignment}.pickle",
    log:
        temp("results/logs/experiment/assigned_counts/assignment/createAssignmentPickleFile.{project}.{assignment}.log"),
    conda:
        getCondaEnv("python3.yaml")
    shell:
        """
        python {input.script} -i {input.files} -o {output} &> {log}
        """


rule experiment_assigned_counts_assignBarcodes:
    """
Assign RNA and DNA barcodes seperately to make the statistic for assigned
"""
    input:
        counts=lambda wc: getFinalCounts(wc.project, wc.config, wc.condition, wc.type, "counts"),
        association="results/experiments/{project}/assignment/{assignment}.tsv.gz",
        script=getScript("count/merge_BC_and_assignment.py"),
    output:
        counts="results/experiments/{project}/assigned_counts/{assignment}/{condition}.{replicate}.{type}.final_counts.config.{config}.tsv.gz",
        statistic=temp(
            "results/experiments/{project}/statistic/assigned_counts/{assignment}/{condition}.{replicate}.{type}.{config}.statistic.tsv.gz"
        ),
    log:
        temp(
            "results/logs/experiment/assigned_counts/assignBarcodes.{project}.{condition}.{replicate}.{type}.{config}.{assignment}.log"
        ),
    conda:
        getCondaEnv("python3.yaml")
    params:
        name="{condition}.{replicate}.{type}",
    shell:
        """
        python {input.script} --counts {input.counts} \
        --assignment {input.association} \
        --output {output.counts} \
        --statistic {output.statistic} \
        --name {params.name} &> {log}
        """


rule experiment_assigned_counts_dna_rna_merge:
    """
Assign merged RNA/DNA barcodes. Filter BC depending on the min_counts option.
"""
    input:
        counts="results/experiments/{project}/counts/{condition}.{replicate}.merged.config.{config}.tsv.gz",
        association="results/experiments/{project}/assignment/{assignment}.tsv.gz",
        script=getScript("count/merge_label.py"),
    output:
        counts="results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}.{replicate}.merged_assigned_counts.tsv.gz",
        bc_counts="results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}.{replicate}.barcode_assigned_counts.tsv.gz",
        removed_bcs="results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}.{replicate}.barcodesRemoved_assigned_counts.tsv.gz",
        statistic=temp(
            "results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}.{replicate}.merged_assigned_counts.statistic.tsv.gz"
        ),
    log:
        temp(
            "results/logs/experiment/assigned_counts/{assignment}/dna_rna_merge.{project}.{condition}.{replicate}.{config}.log"
        ),
    conda:
        getCondaEnv("python3.yaml")
    params:
        minRNACounts=lambda wc: config["experiments"][wc.project]["configs"][wc.config]["filter"]["min_rna_counts"],
        minDNACounts=lambda wc: config["experiments"][wc.project]["configs"][wc.config]["filter"]["min_dna_counts"],
        outlier_detection=lambda wc: (
            "--outlier-detection %s "
            % config["experiments"][wc.project]["configs"][wc.config]["filter"]["outlier_detection"]["method"]
            if "method" in config["experiments"][wc.project]["configs"][wc.config]["filter"]["outlier_detection"]
            else ""
        ),
        outlier_zscore_times=lambda wc: "--outlier-rna-zscore-times %f"
        % config["experiments"][wc.project]["configs"][wc.config]["filter"]["outlier_detection"]["times_zscore"],
    shell:
        """
        python {input.script} --counts {input.counts} \
        --minRNACounts {params.minRNACounts} --minDNACounts {params.minDNACounts} \
        --assignment {input.association} \
        {params.outlier_detection} --outlier-barcodes {output.removed_bcs} \
        {params.outlier_zscore_times} \
        --output {output.counts} \
        --bcOutput {output.bc_counts} \
        --statistic {output.statistic} &> {log}
        """


rule experiment_assigned_counts_make_master_tables:
    """
Final master table with all replicates combined. With and without threshold.
"""
    input:
        counts=lambda wc: expand(
            "results/experiments/{{project}}/assigned_counts/{{assignment}}/{{config}}/{{condition}}.{replicate}.merged_assigned_counts.tsv.gz",
            replicate=getReplicatesOfCondition(wc.project, wc.condition),
        ),
        script=getScript("count/make_master_tables.R"),
    output:
        statistic="results/experiments/{project}/statistic/assigned_counts/{assignment}/{config}/{condition}.average_allreps.merged.tsv.gz",
        all="results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}.allreps.merged.tsv.gz",
        thresh="results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}.allreps_minThreshold.merged.tsv.gz",
    log:
        temp("results/logs/experiment/assigned_counts/make_master_tables.{project}.{condition}.{config}.{assignment}.log"),
    conda:
        getCondaEnv("r.yaml")
    params:
        files=lambda wc: ",".join(
            expand(
                "results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}.{replicate}.merged_assigned_counts.tsv.gz",
                replicate=getReplicatesOfCondition(wc.project, wc.condition),
                project=wc.project,
                condition=wc.condition,
                assignment=wc.assignment,
                config=wc.config,
            )
        ),
        replicates=lambda wc: ",".join(getReplicatesOfCondition(wc.project, wc.condition)),
        thresh=lambda wc: config["experiments"][wc.project]["configs"][wc.config]["filter"]["bc_threshold"],
    shell:
        """
        Rscript {input.script} \
        --threshold {params.thresh} \
        --files {params.files} \
        --replicates {params.replicates} \
        --output {output.thresh} \
        --output-all {output.all} \
        --statistic {output.statistic} &> {log}
        """


rule experiment_assigned_counts_combine_replicates_barcode_output:
    """
Combine replictes of assigned barcode counts into one file.
"""
    input:
        bc_counts=lambda wc: expand(
            "results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}.{replicate}.barcode_assigned_counts.tsv.gz",
            replicate=getReplicatesOfCondition(wc.project, wc.condition),
            project=wc.project,
            condition=wc.condition,
            assignment=wc.assignment,
            config=wc.config,
        ),
        script=getScript("count/merge_replicates_barcode_counts.py"),
    output:
        bc_merged_thresh="results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}.allreps_minThreshold.merged_barcode_assigned_counts.tsv.gz",
        bc_merged_all="results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}.allreps.merged_barcode_assigned_counts.tsv.gz",
    log:
        temp(
            "results/logs/experiment/assigned_counts/combine_replicates_barcode_output.{project}.{condition}.{config}.{assignment}.log"
        ),
    conda:
        getCondaEnv("python3.yaml")
    params:
        thresh=lambda wc: config["experiments"][wc.project]["configs"][wc.config]["filter"]["bc_threshold"],
        bc_counts=lambda wc: " ".join(
            [
                "--counts %s results/experiments/%s/assigned_counts/%s/%s/%s.%s.barcode_assigned_counts.tsv.gz"
                % (rep, wc.project, wc.assignment, wc.config, wc.condition, rep)
                for rep in getReplicatesOfCondition(wc.project, wc.condition)
            ]
        ),
    shell:
        """
        python {input.script} \
        {params.bc_counts} \
        --threshold {params.thresh} \
        --output-threshold {output.bc_merged_thresh} \
        --output {output.bc_merged_all} &> {log}
        """


rule experiment_assigned_counts_combine_replicates:
    """
Combine replicates of master table by summing counts up and using also the average.
"""
    input:
        master_table="results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}.{allreps_or_threshold}.merged.tsv.gz",
        script=getScript("count/combine_replicates.py"),
    output:
        "results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}.{allreps_or_threshold}.merged.combined.tsv.gz",
    log:
        temp(
            "results/logs/experiment/assigned_counts/combine_replicates.{project}.{condition}.{config}.{assignment}.{allreps_or_threshold}.log"
        ),
    conda:
        getCondaEnv("python3.yaml")
    params:
        label_file=lambda wc: (
            "--labels %s" % config["experiments"][wc.project]["label_file"]
            if "label_file" in config["experiments"][wc.project]
            else ""
        ),
    shell:
        """
        python {input.script} \
        --input {input.master_table} \
        {params.label_file} \
        --output {output}  &> {log}
        """


rule experiment_assigned_counts_copy_final_all_files:
    """
Will copy final files to the main folder so that it is clear which files to use.
"""
    input:
        all=lambda wc: "results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}.allreps.merged.tsv.gz",
        bc_all=lambda wc: "results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}.allreps.merged_barcode_assigned_counts.tsv.gz",
    output:
        all="results/experiments/{project}/reporter_experiment.oligo.{condition}.{assignment}.{config}.all.tsv.gz",
        bc_all="results/experiments/{project}/reporter_experiment.barcode.{condition}.{assignment}.{config}.all.tsv.gz",
    log:
        temp("results/logs/experiment/assigned_counts/copy_final_all_files.{project}.{condition}.{assignment}.{config}.log"),
    conda:
        getCondaEnv("default.yaml")
    shell:
        """
        cp {input.all} {output.all} &> {log}
        cp {input.bc_all} {output.bc_all} &>> {log}
        """


rule experiment_assigned_counts_copy_final_thresh_files:
    """
Will copy final files to the main folder so that it is clear which files to use.
"""
    input:
        thresh=lambda wc: "results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}.allreps_minThreshold.merged.tsv.gz",
        bc_thresh=lambda wc: "results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}.allreps_minThreshold.merged_barcode_assigned_counts.tsv.gz",
    output:
        thresh="results/experiments/{project}/reporter_experiment.oligo.{condition}.{assignment}.{config}.min_oligo_threshold_{threshold}.tsv.gz",
        bc_thresh="results/experiments/{project}/reporter_experiment.barcode.{condition}.{assignment}.{config}.min_oligo_threshold_{threshold}.tsv.gz",
    log:
        temp(
            "results/logs/experiment/assigned_counts/copy_final_thresh_files.{project}.{condition}.{assignment}.{config}.{threshold}.log"
        ),
    conda:
        getCondaEnv("default.yaml")
    shell:
        """
        cp {input.thresh} {output.thresh} &> {log}
        cp {input.bc_thresh} {output.bc_thresh} &>> {log}
        """
