.. _Plasmid example:

.. role:: bash(code)
   :language: bash

===========================
ENCODE data (Plasmid based)
===========================

This example runs the assignment and experiment workflow on data from Ryan Tewhey's lab, available via the ENCODE portal. We use the A549 experiment data published in the article `Gosai SJ, Castro RI, Fuentes N et al. Machine-guided design of cell-type-targeting cis-regulatory elements. Nature. 2024. <https://doi.org/10.1038/s41586-024-08070-z>`_. 

The main differences from other examples are:
- For the assignment part, the barcode is attached to the forward read, so we must split them first based on a linker or, in this case, a specific sequence length after the barcode.
- For the experiment/count workflow, as is often done in plasmid-based assays, only one DNA sequencing is performed before electroporation. We need to define the DNA count sequencing fastq file for each replicate. MPRAsnakeflow will recognize that these files are identical, count DNA only once, and use the result for all replicates.

Prerequisites
======================

This example depends on the following data and software:


Installation of MPRAsnakeflow
----------------------------------------

Please install conda, the MPRAsnakeflow environment and clone the actual MPRAsnakeflow master branch. You will find more help under :ref:`Installation`.

Design file
-------------

We need the 

.. code-block:: bash

   mkdir -p ENCFF074MMO/data/fasta
   cd ENCFF074MMO/data/fasta
   wget https://www.encodeproject.org/files/ENCFF074MMO/@@download/ENCFF074MMO.fasta.gz
   cd ../../

Reads association data
-----------------------

There are four sets of association sequencing for this data, which contains a barcode and oligo-forward (FW read) and  oligo-reverse (REV read) sequence. These data must be downloaded. All data is publicly available on ENCODE.

.. note:: You need 23 GB disk space to the download association data!

.. code-block:: bash

    mkdir -p data/oligo_barcode
    cd data/oligo_barcode
    wget https://www.encodeproject.org/files/ENCFF215NWC/@@download/ENCFF215NWC.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF235LUE/@@download/ENCFF235LUE.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF317FNS/@@download/ENCFF317FNS.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF381TUK/@@download/ENCFF381TUK.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF420UFL/@@download/ENCFF420UFL.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF504ZPY/@@download/ENCFF504ZPY.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF708MPZ/@@download/ENCFF708MPZ.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF847HXK/@@download/ENCFF847HXK.fastq.gz
    
    cd ../../

.. note:: Please be sure that all files are downloaded completely without errors! Please check the md5sum of files (see the ENCODE portal).

With

.. code-block:: bash

    tree data


The folder should look like this:

.. code-block:: text

    data
    ├── fasta
    │   └── ENCFF074MMO.fasta.gz
    └── oligo_barcode
        ├── ENCFF215NWC.fastq.gz
        ├── ENCFF235LUE.fastq.gz
        ├── ENCFF317FNS.fastq.gz
        ├── ENCFF381TUK.fastq.gz
        ├── ENCFF420UFL.fastq.gz
        ├── ENCFF504ZPY.fastq.gz
        ├── ENCFF708MPZ.fastq.gz
        └── ENCFF847HXK.fastq.gz


Reads count data
-----------------

Here we have ten sequencing runs for the DNA count data (before electroporation). And then we have five sequencing runs for RNA. For each replicate one run. So five replicates in total.

.. note:: You need 14 GB disk space to the download the count data!

.. code-block:: bash

    mkdir -p data/rna_dna_counts
    cd data/rna_dna_counts
    wget https://www.encodeproject.org/files/ENCFF019RUN/@@download/ENCFF019RUN.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF168OJL/@@download/ENCFF168OJL.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF448RQK/@@download/ENCFF448RQK.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF491IXU/@@download/ENCFF491IXU.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF564JPU/@@download/ENCFF564JPU.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF696HJK/@@download/ENCFF696HJK.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF850RIY/@@download/ENCFF850RIY.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF891CIZ/@@download/ENCFF891CIZ.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF944CEQ/@@download/ENCFF944CEQ.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF966RRE/@@download/ENCFF966RRE.fastq.gz
    

    wget https://www.encodeproject.org/files/ENCFF061UCM/@@download/ENCFF061UCM.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF311KDS/@@download/ENCFF311KDS.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF371LCK/@@download/ENCFF371LCK.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF501AHK/@@download/ENCFF501AHK.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF554OMB/@@download/ENCFF554OMB.fastq.gz

    cd ../../

.. note:: Please be sure that all files are downloaded completely without errors! Please check the md5sum of files (see the ENCODE portal).

With

.. code-block:: bash

    tree data


The folder should look like this:

.. code-block:: text

    data
    ├── fasta
    │   └── ENCFF074MMO.fasta.gz
    ├── oligo_barcode
    │   ├── ENCFF215NWC.fastq.gz
    │   ├── ENCFF235LUE.fastq.gz
    │   ├── ENCFF317FNS.fastq.gz
    │   ├── ENCFF381TUK.fastq.gz
    │   ├── ENCFF420UFL.fastq.gz
    │   ├── ENCFF504ZPY.fastq.gz
    │   ├── ENCFF708MPZ.fastq.gz
    │   └── ENCFF847HXK.fastq.gz
    └── rna_dna_counts
        ├── ENCFF019RUN.fastq.gz
        ├── ENCFF061UCM.fastq.gz
        ├── ENCFF168OJL.fastq.gz
        ├── ENCFF311KDS.fastq.gz
        ├── ENCFF371LCK.fastq.gz
        ├── ENCFF448RQK.fastq.gz
        ├── ENCFF491IXU.fastq.gz
        ├── ENCFF501AHK.fastq.gz
        ├── ENCFF554OMB.fastq.gz
        ├── ENCFF564JPU.fastq.gz
        ├── ENCFF696HJK.fastq.gz
        ├── ENCFF850RIY.fastq.gz
        ├── ENCFF891CIZ.fastq.gz
        ├── ENCFF944CEQ.fastq.gz
        └── ENCFF966RRE.fastq.gz




MPRAsnakeflow
=================================

We will run assignmenta nd count workflow together. But it is of course possible to run them seperately using different config files. Then you have to use the assignment `fromFile` not `fromConfig`. But first we need to define the config file and the experiment CSV file to map DNA/RNA counts to the correct replicate.


Create config files
-------------------

.. code-block:: bash

    cat << 'EOF' >  config.yaml
    ---
    version: "0.5.4"
    assignments:
        ENCFF074MMOAssignment:
            bc_length: 20
            BC_rev_comp: false
            linker: TCTAGAGGTTCGTCGACGCGATCGCAGGAGCCGCAGTG
            adapters:
                3prime:
                    - CGTCAAGCGGCCAGTT
            alignment_tool:
                split_number: 30
                tool: bbmap
                configs:
                    sequence_length: 200
                    alignment_start: 1
            FW:
                - data/oligo_barcode/ENCFF235LUE.fastq.gz
                - data/oligo_barcode/ENCFF708MPZ.fastq.gz
                - data/oligo_barcode/ENCFF381TUK.fastq.gz
                - data/oligo_barcode/ENCFF847HXK.fastq.gz
            REV:
                - data/oligo_barcode/ENCFF215NWC.fastq.gz
                - data/oligo_barcode/ENCFF317FNS.fastq.gz
                - data/oligo_barcode/ENCFF504ZPY.fastq.gz
                - data/oligo_barcode/ENCFF420UFL.fastq.gz
            design_file: data/fasta/ENCFF074MMO.fasta.gz
            configs:
                default: {}
    experiments:
        ENCFF074MMOExperiment:
            bc_length: 20
            data_folder: data/rna_dna_counts
            experiment_file: experiment.csv
            assignments:
                ENCFF074MMOAssignment:
                    type: config
                    assignment_name: ENCFF074MMOAssignment
                    assignment_config: default
            configs:
                default: {}
    EOF

And the :code:`experiment.csv` file to map the DNA/RNA counts to the correct replicate. The experiment file is a simple CSV file with the following content:

.. code-block:: bash

    cat << 'EOF' >  experiment.csv
    Condition,Replicate,DNA_BC_F,RNA_BC_F
    A549,1,ENCFF448RQK.fastq.gz;ENCFF019RUN.fastq.gz;ENCFF850RIY.fastq.gz;ENCFF966RRE.fastq.gz;ENCFF168OJL.fastq.gz;ENCFF891CIZ.fastq.gz;ENCFF696HJK.fastq.gz;ENCFF491IXU.fastq.gz;ENCFF944CEQ.fastq.gz;ENCFF564JPU.fastq.gz,ENCFF311KDS.fastq.gz
    A549,2,ENCFF448RQK.fastq.gz;ENCFF019RUN.fastq.gz;ENCFF850RIY.fastq.gz;ENCFF966RRE.fastq.gz;ENCFF168OJL.fastq.gz;ENCFF891CIZ.fastq.gz;ENCFF696HJK.fastq.gz;ENCFF491IXU.fastq.gz;ENCFF944CEQ.fastq.gz;ENCFF564JPU.fastq.gz,ENCFF554OMB.fastq.gz
    A549,3,ENCFF448RQK.fastq.gz;ENCFF019RUN.fastq.gz;ENCFF850RIY.fastq.gz;ENCFF966RRE.fastq.gz;ENCFF168OJL.fastq.gz;ENCFF891CIZ.fastq.gz;ENCFF696HJK.fastq.gz;ENCFF491IXU.fastq.gz;ENCFF944CEQ.fastq.gz;ENCFF564JPU.fastq.gz,ENCFF501AHK.fastq.gz
    A549,4,ENCFF448RQK.fastq.gz;ENCFF019RUN.fastq.gz;ENCFF850RIY.fastq.gz;ENCFF966RRE.fastq.gz;ENCFF168OJL.fastq.gz;ENCFF891CIZ.fastq.gz;ENCFF696HJK.fastq.gz;ENCFF491IXU.fastq.gz;ENCFF944CEQ.fastq.gz;ENCFF564JPU.fastq.gz,ENCFF371LCK.fastq.gz
    A549,5,ENCFF448RQK.fastq.gz;ENCFF019RUN.fastq.gz;ENCFF850RIY.fastq.gz;ENCFF966RRE.fastq.gz;ENCFF168OJL.fastq.gz;ENCFF891CIZ.fastq.gz;ENCFF696HJK.fastq.gz;ENCFF491IXU.fastq.gz;ENCFF944CEQ.fastq.gz;ENCFF564JPU.fastq.gz,ENCFF061UCM.fastq.gz
        EOF

Run snakemake
-------------

Now we are ready to run MPRAsnakeflow. We will do it on one node with 50GB memory and 30 cores.

We will run the pipeline directly in the :code:`ENCFF074MMO` folder. The MPRAsnakeflow workflow can be in a different directory. Let's assume :code:`/home/user/MPRAsnakeflow`. 

First we do a try run using snakemake :code:`-n` option. The MPRAsnakeflow command is:

.. code-block:: bash

    conda activate mprasnakeflow
    snakemake -c 1 --sdm apptainer conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile config.yaml -n --quiet rules

You should see a list of rules that will be executed. Here is the summary:

.. code-block:: text

    Job stats:                                                                      
    job                                                                        count
    -----------------------------------------------------------------------  -------
    all  1
    assignment_3prime_remove                                                      30
    assignment_attach_idx                                                         60
    assignment_check_design                                                        1
    assignment_collect                                                             1
    assignment_collectBCs                                                          1
    assignment_fastq_split                                                         3
    assignment_filter                                                              1
    assignment_flagstat                                                            1
    assignment_hybridFWRead_get_reads_by_cutadapt                                  1
    assignment_idx_bam                                                             1
    assignment_mapping_bbmap                                                      30
    assignment_mapping_bbmap_getBCs                                               30
    assignment_merge                                                              30
    assignment_statistic_assignedCounts                                            1
    assignment_statistic_assignment                                                1
    assignment_statistic_quality_metric                                            1
    assignment_statistic_totalCounts                                               1
    experiment_assigned_counts_assignBarcodes                                     10
    experiment_assigned_counts_combine_replicates                                  2
    experiment_assigned_counts_combine_replicates_barcode_output                   1
    experiment_assigned_counts_copy_final_all_files                                1
    experiment_assigned_counts_copy_final_thresh_files                             1
    experiment_assigned_counts_dna_rna_merge                                       5
    experiment_assigned_counts_filterAssignment                                    1
    experiment_assigned_counts_make_master_tables                                  1
    experiment_counts_dna_rna_merge_counts                                        10
    experiment_counts_filter_counts                                                6
    experiment_counts_final_counts                                                 6
    experiment_counts_onlyFW_raw_counts_by_length                                  6
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
    experiment_statistic_counts_BC_in_RNA_DNA                                     10
    experiment_statistic_counts_BC_in_RNA_DNA_merge                                2
    experiment_statistic_counts_barcode_base_composition                           6
    experiment_statistic_counts_final                                              2
    experiment_statistic_counts_frequent_umis                                      6
    experiment_statistic_counts_stats_merge                                        2
    experiment_statistic_counts_table                                             12
    experiment_statistic_quality_metric                                            1
    qc_report_assoc                                                                1
    qc_report_count                                                                1
    total                                                                        307


When dry-run does not give any errors we will run the workflow. We use a machine with 30 threads/cores to run the workflow and 60GB memory. Therefore :code:`split_number` is set to 30 to parallelize the workflow. Also we are using 10 threads for mapping (bbmap). But snakemake takes care that no more than 30 threads are used.

.. code-block:: bash

    snakemake -c 30 --sdm conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile config.yaml -q --set-threads assignment_mapping_bbmap=10 --resources mem_mb=60000


.. note:: Please modify your code when running in a cluster environment. We have an example SLURM profile within the MPRAsnakeflow repository under :code:`profiles/default/config.yaml`. You can use it in snakemake with :code:`--workflow-profile $PIPELINE/profiles/default`. But adapt your before the :code:`slurm_partition`

Results
-----------------

For the assignment all output files will be in the :code:`results/assignment/ENCFF074MMOAssignment` folder. The final assignment is in :code:`results/assignment/ENCFF074MMOAssignment/assignment_barcodes.default.tsv.gz`. Also you should have a look at the qc report: :code:`results/assignment/ENCFF074MMOAssignment/qc_report.default.html`. You can find an example qc report here: `Example assignment QC report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/plasmid_assignment.qc_report.default.html>`_.


For the experiment all output files will be in the :code:`results/experiment/ENCFF074MMOExperiment` folder. The final count files is :code:`results/experiment/ENCFF074MMOExperiment/reporter_experiment.barcode.A549.ENCFF074MMOAssignment.default.all.tsv.gz` for the barcode file and :code:`results/experiment/ENCFF074MMOExperiment/reporter_experiment.oligo.A549.ENCFF074MMOAssignment.default.all.tsv.gz` for the aggregated oligo files. Also you should have a look at the qc report: :code:`results/experiment/ENCFF074MMOExperiment/qc_report.A549.ENCFF074MMOAssignment.default.html`. You can find an example qc report here: `Example experiment QC report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/plasmid_experiment.qc_report.A549.ENCFF074MMOAssignment.default.html>`_.
