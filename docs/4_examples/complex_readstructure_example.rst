.. _Complex example:

.. role:: bash(code)
    :language: bash

======================
Complex read structure
======================

This example runs the assignment and experiment workflow on data from the publication `Nathan S. Abell et al. Multiple causal variants underlie genetic associations in humans. Science 375, 1247-1254 (2022). <https://doi.org/10.1126/science.abj5117>`_. Their data was published on `GEO:GSE174534 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174534>`_. They used a custom processing pipeline, available on `Zenodo <https://zenodo.org/records/5921042>`_, but we try to run everything in MPRAsnakeflow.

Looking at the data, we see that the reads contain more information than only the barcodes or barcode and oligo information that is needed for the pipeline. So we have to enable the MPRAsnakeflow pre-processing steps. Unfortunately the two runs of the uploaded assignment data are already different (151 bp paired-end reads for SRR14567164 and 150 bp for SRR14567165). Because MPRAsnakeflow assumes that all input reads are of the same structure, we have to trim the first run by one base. In this case it is only necessary for the FWD read. The REV read can stay with 151 bp because the oligo overlap is just 1 bp longer.

Another solvable issue is that the oligo sequences are identical (we simply take each sequence only once) and oligos are in reverse complement orientation in the reads. This is no problem for the mapper (we map in both directions), but because sequences with forward and reverse orientation are in the dataset, we have to enable the strand sensitive option via:

.. code-block:: yaml

     strand_sensitive:
        enable: true

Within the strand sensitive mode, reads get a unique sequence at the end of the read as well as the reference oligo file. Unfortunately, because oligo sequences are sequenced in reverse complement, this does not work. But a simple workaround is to use the input oligo design as reverse complement.


.. note:: This example is very large and you should have around 800GB disk space available to download and process the data!


Prerequisites
=============

This example depends on the following data and software:


Installation of MPRAsnakeflow
------------------------------

Please install conda, the MPRAsnakeflow environment and clone the actual MPRAsnakeflow master branch. You will find more help under :ref:`Installation`. In addition install `cutadapt`, `seqkit` and `sra-toolkit` e.g. via conda:

.. code-block:: bash

    conda install -c bioconda cutadapt seqkit sra-tools

Design file
-----------

We need the design file and have to modify it by trimming not sequenced adapters as well as creating the reverse complement (see above). 

.. code-block:: bash

    mkdir -p data/association
    cd data/association
    wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5319012&format=file&file=GSM5319012%5FVariant%2DOligo%2DDesign%2Etxt%2Egz" -O GSM5319012_Variant-Oligo-Design.txt.gz
    # unique sequence and trimming
    zcat GSM5319012_Variant-Oligo-Design.txt.gz | awk '{print substr($2, 16, length($2)-30)}' | sort | uniq | awk '{print ">ID"NR"\n"$1}' > variant_oligo_design_trimmed.fa
    # reverse complement
    seqkit seq --seq-type DNA -w 200 --reverse --complement variant_oligo_design_trimmed.fa > variant_oligo_design_trimmed_RC.fa
    cd ../../

Reads association data
----------------------

There are four sets of association sequencing for this data, which contain a barcode and oligo-forward (FW read) and oligo-reverse (REV read) sequence. These data must be downloaded. All data is publicly available on ENCODE.

.. code-block:: bash

    mkdir -p data/association
    cd data/association
    prefetch --max-size 50GB SRR14567165
    fastq-dump --gzip --split-files SRR14567165

    prefetch --max-size 200GB SRR14567164
    fastq-dump --gzip --split-files SRR14567164

    # trim 1 bp from FWD read of SRR14567164
    cutadapt -u 1 <(zcat SRR14567164_1.fastq.gz) | gzip -c > SRR14567164_1.trim1bp.fastq.gz
        
    cd ../../

.. note:: Please be sure that all files are downloaded completely without errors! Please check the md5sum of files (see the ENCODE portal).

With

.. code-block:: bash

    tree data


The folder should look like this:

.. code-block:: text

    data/
    └── association
        ├── GSM5319012_Variant-Oligo-Design.txt.gz
        ├── SRR14567164_1.fastq.gz
        ├── SRR14567164_1.trim1bp.fastq.gz
        ├── SRR14567164_2.fastq.gz
        ├── SRR14567165_1.fastq.gz
        ├── SRR14567165_2.fastq.gz
        ├── variant_oligo_design_trimmed.fa
        └── variant_oligo_design_trimmed_RC.fa

Reads count data
----------------

We have three replicates of DNA and RNA counts each. These data must be downloaded.

replicate | plasmid | cDNA 
1   | SRR14567158   | SRR14567161
2   | SRR14567159   | SRR14567162
3   | SRR14567160   | SRR14567163

.. code-block:: bash

    mkdir -p data/counts
    cd data/counts
    # download data
    prefetch --max-size 50GB SRR14567158 SRR14567159 SRR14567160 SRR14567161 SRR14567162 SRR14567163
    fastq-dump --gzip --split-files SRR14567158 SRR14567159 SRR14567160 SRR14567161 SRR14567162 SRR14567163

    cd ../../

With

.. code-block:: bash

    tree data


The folder should look like this:

.. code-block:: text

    data/
    ├── association
    │   ├── GSM5319012_Variant-Oligo-Design.txt.gz
    │   ├── SRR14567164_1.fastq.gz
    │   ├── SRR14567164_1.trim1bp.fastq.gz
    │   ├── SRR14567164_2.fastq.gz
    │   ├── SRR14567165_1.fastq.gz
    │   ├── SRR14567165_2.fastq.gz
    │   ├── variant_oligo_design_trimmed.fa
    │   └── variant_oligo_design_trimmed_RC.fa
    └── counts
        ├── SRR14567158_1.fastq.gz
        ├── SRR14567159_1.fastq.gz
        ├── SRR14567160_1.fastq.gz
        ├── SRR14567161_1.fastq.gz
        ├── SRR14567162_1.fastq.gz
        └── SRR14567163_1.fastq.gz




MPRAsnakeflow
=============

We will run assignment and count workflow together. But it is of course possible to run them separately using different config files. Then you have to use the assignment `fromFile` not `fromConfig`. But first we need to define the config file and the experiment CSV file to map DNA/RNA counts to the correct replicate.

Some important information you have to know when creating the config file. The oligo design is 150bp and the barcodes are 20bp long. We have to remove a 16bp sequence (seems to be random) from the beginning of the reverse read (or the 3' end later of a merged read). The linker sequence between barcode and oligo is 38bp long. We have to enable strand sensitive mode.

For the counts the reads are longer than the 20bp barcode. But for single-end reads we always use the first N bases (depending on the BC length) as a barcode. Therefore no adapter removal is necessary for the counts.

Create config files
-------------------

.. code-block:: bash

    cat << 'EOF' >  config.yaml
    ---
    version: "0.6"
    assignments:
        SRR14567164SRR14567165:
            bc_length: 20
            linker_length: 38
            adapters:
                REV:
                    - 16
            strand_sensitive:
                enable: true
            alignment_tool:
                split_number: 30
                tool: bbmap
                configs:
                    sequence_length: 150
                    alignment_start: 1
            FWD:
                - data/association/SRR14567164_1.trim1bp.fastq.gz
                - data/association/SRR14567165_1.fastq.gz
            REV:
                - data/association/SRR14567164_2.fastq.gz
                - data/association/SRR14567165_2.fastq.gz
            design_file: data/association/variant_oligo_design_trimmed_RC.fa
            configs:
                default: {}
    experiments:
        abellEtAl:
            bc_length: 20
            data_folder: data/counts
            experiment_file: experiment.csv
            assignments:
                SRR14567164SRR14567165Assignment:
                    type: config
                    assignment_name: SRR14567164SRR14567165
                    assignment_config: default
            configs:
                default: {}
    EOF

And the :code:`experiment.csv` file to map the DNA/RNA counts to the correct replicate. The experiment file is a simple CSV file with the following content:

.. code-block:: bash

    cat << 'EOF' >  experiment.csv
    Condition,Replicate,DNA_BC_F,RNA_BC_F
    GM12878,1,SRR14567158_1.fastq.gz,SRR14567161_1.fastq.gz
    GM12878,2,SRR14567159_1.fastq.gz,SRR14567162_1.fastq.gz
    GM12878,3,SRR14567160_1.fastq.gz,SRR14567163_1.fastq.gz
    EOF


Run snakemake
-------------

Now we are ready to run MPRAsnakeflow. We will do it on one node with 50GB memory and 30 cores.

We will run the pipeline directly in the actual folder. The MPRAsnakeflow workflow can be in a different directory. Let's assume :code:`/home/user/MPRAsnakeflow`. 

First we do a dry run using snakemake :code:`-n` option. The MPRAsnakeflow command is:

.. code-block:: bash

    conda activate mprasnakeflow
    snakemake -c 1 --sdm apptainer conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile config.yaml -n --quiet rules

You should see a list of rules that will be executed. Here is the summary:

.. code-block:: text

    Job stats:
    job                                                                        count
    -----------------------------------------------------------------------  -------
    all                                                                            1
    assignment_attach_idx                                                         60
    assignment_check_design                                                        1
    assignment_collect                                                             1
    assignment_collectBCs                                                          1
    assignment_fastq_split                                                         3
    assignment_filter                                                              1
    assignment_flagstat                                                            1
    assignment_hybridFWDRead_get_reads_by_length                                   1
    assignment_idx_bam                                                             1
    assignment_mapping_bbmap                                                      30
    assignment_mapping_bbmap_getBCs                                               30
    assignment_merge                                                              30
    assignment_preprocessingadapter_remove                                         1
    assignment_statistic_assignedCounts                                            1
    assignment_statistic_assignment                                                1
    assignment_statistic_quality_metric                                            1
    assignment_statistic_totalCounts                                               1
    experiment_assigned_counts_assignBarcodes                                      6
    experiment_assigned_counts_combine_replicates                                  2
    experiment_assigned_counts_combine_replicates_barcode_output                   1
    experiment_assigned_counts_copy_final_all_files                                1
    experiment_assigned_counts_copy_final_thresh_files                             1
    experiment_assigned_counts_dna_rna_merge                                       3
    experiment_assigned_counts_filterAssignment                                    1
    experiment_assigned_counts_make_master_tables                                  1
    experiment_counts_dna_rna_merge_counts                                         6
    experiment_counts_filter_counts                                                6
    experiment_counts_final_counts                                                 6
    experiment_counts_onlyFWD_raw_counts                                           6
    experiment_statistic_assigned_counts_combine_BC_assignment_stats               1
    experiment_statistic_assigned_counts_combine_BC_assignment_stats_helper        1
    experiment_statistic_assigned_counts_combine_stats_dna_rna_merge               1
    experiment_statistic_assigned_counts_combine_stats_dna_rna_merge_all           1
    experiment_statistic_bc_overlap_combine_assigned_counts                        1
    experiment_statistic_bc_overlap_combine_counts                                 1
    experiment_statistic_bc_overlap_run                                            4
    experiment_statistic_correlation_bc_counts                                     2
    experiment_statistic_correlation_bc_counts_hist                                2
    experiment_statistic_correlation_calculate                                     1
    experiment_statistic_correlation_combine_bc_assigned                           1
    experiment_statistic_correlation_combine_bc_raw                                1
    experiment_statistic_correlation_combine_oligo                                 1
    experiment_statistic_correlation_hist_box_plots                                1
    experiment_statistic_counts_BC_in_RNA_DNA                                      6
    experiment_statistic_counts_BC_in_RNA_DNA_merge                                2
    experiment_statistic_counts_barcode_base_composition                           6
    experiment_statistic_counts_final                                              2
    experiment_statistic_counts_frequent_umis                                      6
    experiment_statistic_counts_stats_merge                                        2
    experiment_statistic_counts_table                                             12
    experiment_statistic_quality_metric                                            1
    qc_report_assoc                                                                1
    qc_report_count                                                                1
    total                                                                        264


When dry-run does not give any errors we will run the workflow. We use a machine with 30 threads/cores to run the workflow and 60GB memory. Therefore :code:`split_number` is set to 30 to parallelize the workflow. Also we are using 10 threads for mapping (bbmap). But snakemake takes care that no more than 30 threads are used.

.. code-block:: bash

    snakemake -c 30 --sdm conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile config.yaml -q --set-threads assignment_mapping_bbmap=10 --resources mem_mb=60000


.. note:: Please modify your code when running in a cluster environment. We have an example SLURM profile within the MPRAsnakeflow repository under :code:`profiles/default/config.yaml`. You can use it in snakemake with :code:`--workflow-profile $PIPELINE/profiles/default`. But adapt it before the :code:`slurm_partition`

Results
-------

For the assignment all output files will be in the :code:`results/assignment/SRR14567164SRR14567165` folder. The final assignment is in :code:`results/assignment/SRR14567164SRR14567165/assignment_barcodes.default.tsv.gz`. Also you should have a look at the QC report: :code:`results/assignment/SRR14567164SRR14567165/qc_report.default.html`. You can find an example QC report here: `Example assignment QC report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/SRR14567164SRR14567165.qc_report.default.html>`_.


For the experiment all output files will be in the :code:`results/experiments/abellEtAl` folder. The final count file is :code:`results/experiments/abellEtAl/reporter_experiment.barcode.GM12878.SRR14567164SRR14567165Assignment.default.all.tsv.gz` for the barcode file and :code:`results/experiments/abellEtAl/reporter_experiment.oligo.GM12878.SRR14567164SRR14567165Assignment.default.all.tsv.gz` for the aggregated oligo files. Also you should have a look at the QC report: :code:`results/experiments/abellEtAl/qc_report.GM12878.SRR14567164SRR14567165Assignment.default.html`. You can find an example QC report here: `Example experiment QC report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/abeletall.qc_report.GM12878.SRR14567164SRR14567165Assignment.default.html>`_.
