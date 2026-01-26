.. _Combined example:

.. role:: bash(code)
   :language: bash

============================
Combined Workflow
============================

This example runs the assignment and the experiment/count workflow on 5'/5' WT MPRA data in the HepG2 cell line from `Klein J., Agarwal, V., Keith, A., et al. 2019 <https://www.biorxiv.org/content/10.1101/576405v1.full.pdf>`_.

Prerequisites
=============

This example depends on the following data and software:

Installation of MPRAsnakeflow
-----------------------------
Please install conda, the MPRAsnakeflow environment, and clone the actual MPRAsnakeflow master branch. You will find more help under :ref:`Installation`.

Meta Data
---------
It is necessary to get the ordered oligo array so that each enhancer sequence can be labeled in the analysis and to trim any adaptors still in the sequence. In this case, we trim off 15bp from the end of each sequence.

.. code-block:: bash

    mkdir -p combined_basic/data
    cd combined_basic/data
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4237nnn/GSM4237954/suppl/GSM4237954_9MPRA_elements.fa.gz

    zcat GSM4237954_9MPRA_elements.fa.gz | awk '{ count+=1; if (count == 1) { print } else { print substr($1,1,171)}; if (count == 2) { count=0 } }' | awk '{gsub(/[\]\[]/,"_")} $0' > design.fa

Reads
-----
There is one set of association sequencing for this data, which contains a forward (CRS-forward), reverse (CRS-reverse), and index (barcode) read for DNA and RNA. These data must be downloaded. All data is publicly available on the Short Read Archive (SRA). We will use SRA-toolkit to obtain the data.

.. note:: You need 10 GB of disk space to download the data!

.. code-block:: bash

    conda install sra-tools
    cd combined_basic/data
    fastq-dump --gzip --split-files SRR10800986
    cd ..

For large files and unstable internet connections, we recommend the command `prefetch` from SRA tools before running `fastq-dump`. This command provides better warnings when something goes wrong.

.. code-block:: bash

    conda install sra-tools
    cd combined_basic/data
    prefetch SRR10800986
    fastq-dump --gzip --split-files SRR10800986
    cd ..

.. note:: Please ensure that all files are downloaded completely without errors! Depending on your internet connection, this can take a while. If you just want some data to run MPRAsnakeflow, you can limit yourself to one condition and/or just one replicate.

With the following command:

.. code-block:: bash

    tree data

The folder should look like this:

.. code-block:: text

    data
    ├── design.fa
    ├── SRR10800986_1.fastq.gz
    ├── SRR10800986_2.fastq.gz
    └── SRR10800986_3.fastq.gz

MPRAsnakeflow
=============
Now we are ready to run MPRAsnakeflow and create CRS-barcode mappings and counts.

Run Snakemake
-------------
Now we have everything at hand to run the combined MPRAsnakeflow pipeline. We will run the pipeline directly in the :code:`combined_basic` folder. The MPRAsnakeflow workflow can be in a different directory. Let's assume :code:`/home/user/MPRAsnakeflow`.

First, configure the config file:

.. literalinclude:: ../../resources/combined_basic/config.yml
   :language: yaml

Perform a dry-run using Snakemake's :code:`-n` option. The MPRAsnakeflow command is:

.. code-block:: bash

    cd combined_basic
    conda activate mprasnakeflow
    snakemake -c 30 --use-conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile /home/user/MPRAsnakeflow/resources/combined_basic/config.yml -n -q --set-threads assignment_mapping_bbmap=10  --resources mem_mb=60000

You should see a list of rules that will be executed. This is the summary:

.. code-block:: text

    Job stats:                                                                                                                                                                                                     job                                                                        count
    -----------------------------------------------------------------------  -------
    all                                                                            1
    assignment_attach_idx                                                         60
    assignment_check_design                                                        1
    assignment_collect                                                             1
    assignment_collectBCs                                                          1
    assignment_fastq_split                                                         3
    assignment_filter                                                              1
    assignment_flagstat                                                            1
    assignment_idx_bam                                                             1
    assignment_mapping_bbmap                                                      30
    assignment_mapping_bbmap_getBCs                                               30
    assignment_merge                                                              30
    experiment_assigned_counts_combine_replicates_barcode_output                   1
    experiment_assigned_counts_copy_final_all_files                                1
    experiment_assigned_counts_copy_final_thresh_files                             1
    experiment_assigned_counts_dna_rna_merge                                       3
    experiment_assigned_counts_filterAssignment                                    1
    experiment_assigned_counts_make_master_tables                                  1
    experiment_counts_dna_rna_merge_counts                                         6
    experiment_counts_filter_counts                                                6
    experiment_counts_final_counts                                                 6
    experiment_counts_umi_create_BAM                                               6
    experiment_counts_umi_raw_counts                                               6
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
    total                                                                        268



When the dry-run does not give any errors, run the workflow. We use a machine with 30 threads/cores to run the workflow. Therefore, :code:`split_number` is set to 30 to parallelize the workflow. Also, we are using 10 threads for mapping (BBMap). Snakemake ensures that no more than 30 threads are used.

.. code-block:: bash

    snakemake -c 30 --use-conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile /home/user/MPRAsnakeflow/resources/combined_basic/config.yml --set-threads assignment_mapping_bbmap=10  --resources mem_mb=60000

Results
-------

All needed output files will be in the :code:`results/assignment/assocBasic` folder for assignment results. The folder :code:`results/experiments/countBasic` contains the count results. A nice overview (QC report) is shown in :code:`results/experiments/countBasic/qc_report.default.html`. This HTML report contains information about statistics tables and plots. You can find an example qc report here: `Example assignment QC report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/combined_example1.assignment.qc_report.default.html>`_ and here `Example Experiment QC report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/combined_example1.count.qc_report.default.html>`_.

To generate a final report, use the following command:

.. code-block:: bash

    snakemake --config config.yml --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --report report.html

This HTML report contains information about the Snakemake run and integrates statistics tables and plots.

Total file tree of the results folder:

.. code-block:: text

    results/
    ├── assignment
    │   └── assocBasic
    │       ├── aligned_merged_reads.bam
    │       ├── aligned_merged_reads.bam.bai
    │       ├── assignment_barcodes.default.tsv.gz
    │       ├── assignment_barcodes_with_ambiguous.default.tsv.gz
    │       ├── barcodes_incl_other.tsv.gz
    │       ├── design_check.done
    │       ├── design_check.err
    │       ├── qc_metrics.default.json
    │       ├── qc_report.default.html
    │       ├── reference
    │       │   └── reference.fa
    │       └── statistic
    │           ├── assigned_counts.default.tsv
    │           ├── assignment
    │           │   └── bam_stats.txt
    │           ├── assignment.default.png
    │           ├── assignment.default.tsv.gz
    │           └── total_counts.tsv
    ├── experiments
    │   └── countBasic
    │       ├── assigned_counts
    │       │   └── fromWorkflow
    │       │       ├── HEPG2.1.DNA.final_counts.config.default.tsv.gz
    │       │       ├── HEPG2.1.RNA.final_counts.config.default.tsv.gz
    │       │       ├── HEPG2.1.merged.config.default.tsv.gz
    │       │       ├── HEPG2.2.DNA.final_counts.config.default.tsv.gz
    │       │       ├── HEPG2.2.RNA.final_counts.config.default.tsv.gz
    │       │       ├── HEPG2.2.merged.config.default.tsv.gz
    │       │       ├── HEPG2.3.DNA.final_counts.config.default.tsv.gz
    │       │       ├── HEPG2.3.RNA.final_counts.config.default.tsv.gz
    │       │       ├── HEPG2.3.merged.config.default.tsv.gz
    │       │       └── default
    │       │           ├── HEPG2.1.barcode_assigned_counts.tsv.gz
    │       │           ├── HEPG2.1.barcodesRemoved_assigned_counts.tsv.gz
    │       │           ├── HEPG2.1.merged_assigned_counts.tsv.gz
    │       │           ├── HEPG2.2.barcode_assigned_counts.tsv.gz
    │       │           ├── HEPG2.2.barcodesRemoved_assigned_counts.tsv.gz
    │       │           ├── HEPG2.2.merged_assigned_counts.tsv.gz
    │       │           ├── HEPG2.3.barcode_assigned_counts.tsv.gz
    │       │           ├── HEPG2.3.barcodesRemoved_assigned_counts.tsv.gz
    │       │           ├── HEPG2.3.merged_assigned_counts.tsv.gz
    │       │           ├── HEPG2.allreps.merged.combined.tsv.gz
    │       │           ├── HEPG2.allreps.merged.tsv.gz
    │       │           ├── HEPG2.allreps.merged_barcode_assigned_counts.tsv.gz
    │       │           ├── HEPG2.allreps_minThreshold.merged.combined.tsv.gz
    │       │           ├── HEPG2.allreps_minThreshold.merged.tsv.gz
    │       │           └── HEPG2.allreps_minThreshold.merged_barcode_assigned_counts.tsv.gz
    │       ├── assignment
    │       │   └── fromWorkflow.tsv.gz
    │       ├── counts
    │       │   ├── HEPG2.1.DNA.filtered_counts.tsv.gz
    │       │   ├── HEPG2.1.DNA.final_counts.tsv.gz
    │       │   ├── HEPG2.1.RNA.filtered_counts.tsv.gz
    │       │   ├── HEPG2.1.RNA.final_counts.tsv.gz
    │       │   ├── HEPG2.1.merged.config.default.tsv.gz
    │       │   ├── HEPG2.2.DNA.filtered_counts.tsv.gz
    │       │   ├── HEPG2.2.DNA.final_counts.tsv.gz
    │       │   ├── HEPG2.2.RNA.filtered_counts.tsv.gz
    │       │   ├── HEPG2.2.RNA.final_counts.tsv.gz
    │       │   ├── HEPG2.2.merged.config.default.tsv.gz
    │       │   ├── HEPG2.3.DNA.filtered_counts.tsv.gz
    │       │   ├── HEPG2.3.DNA.final_counts.tsv.gz
    │       │   ├── HEPG2.3.RNA.filtered_counts.tsv.gz
    │       │   ├── HEPG2.3.RNA.final_counts.tsv.gz
    │       │   ├── HEPG2.3.merged.config.default.tsv.gz
    │       │   ├── useUMI.HEPG2.1.DNA.bam
    │       │   ├── useUMI.HEPG2.1.DNA.raw_counts.tsv.gz
    │       │   ├── useUMI.HEPG2.1.RNA.bam
    │       │   ├── useUMI.HEPG2.1.RNA.raw_counts.tsv.gz
    │       │   ├── useUMI.HEPG2.2.DNA.bam
    │       │   ├── useUMI.HEPG2.2.DNA.raw_counts.tsv.gz
    │       │   ├── useUMI.HEPG2.2.RNA.bam
    │       │   ├── useUMI.HEPG2.2.RNA.raw_counts.tsv.gz
    │       │   ├── useUMI.HEPG2.3.DNA.bam
    │       │   ├── useUMI.HEPG2.3.DNA.raw_counts.tsv.gz
    │       │   ├── useUMI.HEPG2.3.RNA.bam
    │       │   └── useUMI.HEPG2.3.RNA.raw_counts.tsv.gz
    │       ├── qc_metrics.HEPG2.fromWorkflow.default.json
    │       ├── qc_report.HEPG2.fromWorkflow.default.html
    │       ├── reporter_experiment.barcode.HEPG2.fromWorkflow.default.all.tsv.gz
    │       ├── reporter_experiment.barcode.HEPG2.fromWorkflow.default.min_oligo_threshold_10.tsv.gz
    │       ├── reporter_experiment.oligo.HEPG2.fromWorkflow.default.all.tsv.gz
    │       ├── reporter_experiment.oligo.HEPG2.fromWorkflow.default.min_oligo_threshold_10.tsv.gz
    │       └── statistic
    │           ├── assigned_counts
    │           │   └── fromWorkflow
    │           │       └── default
    │           │           ├── HEPG2.DNA.pairwise.minThreshold.png
    │           │           ├── HEPG2.DNA.pairwise.png
    │           │           ├── HEPG2.RNA.pairwise.minThreshold.png
    │           │           ├── HEPG2.RNA.pairwise.png
    │           │           ├── HEPG2.Ratio.pairwise.minThreshold.png
    │           │           ├── HEPG2.Ratio.pairwise.png
    │           │           ├── HEPG2.average_allreps.merged.tsv.gz
    │           │           ├── HEPG2.barcodesPerInsert.png
    │           │           ├── HEPG2.dna_vs_rna.png
    │           │           ├── HEPG2.dna_vs_rna_minThreshold.png
    │           │           ├── HEPG2.group_barcodesPerInsert_box.png
    │           │           └── HEPG2.group_barcodesPerInsert_box_minThreshold.png
    │           ├── barcode
    │           │   ├── assigned_counts
    │           │   │   └── fromWorkflow
    │           │   │       ├── HEPG2.default.DNA.perBarcode.png
    │           │   │       ├── HEPG2.default.RNA.perBarcode.png
    │           │   │       ├── HEPG2.default.barcode.DNA.pairwise.png
    │           │   │       ├── HEPG2.default.barcode.RNA.pairwise.png
    │           │   │       └── HEPG2.default.barcode.Ratio.pairwise.png
    │           │   └── counts
    │           │       ├── HEPG2.default.DNA.perBarcode.png
    │           │       ├── HEPG2.default.RNA.perBarcode.png
    │           │       ├── HEPG2.default.barcode.DNA.pairwise.png
    │           │       ├── HEPG2.default.barcode.RNA.pairwise.png
    │           │       └── HEPG2.default.barcode.Ratio.pairwise.png
    │           ├── bc_overlap.assigned_counts.default.fromWorkflow.tsv
    │           ├── bc_overlap.counts.default.tsv
    │           ├── counts
    │           │   ├── BCNucleotideComposition.HEPG2.1.DNA.tsv.gz
    │           │   ├── BCNucleotideComposition.HEPG2.1.RNA.tsv.gz
    │           │   ├── BCNucleotideComposition.HEPG2.2.DNA.tsv.gz
    │           │   ├── BCNucleotideComposition.HEPG2.2.RNA.tsv.gz
    │           │   ├── BCNucleotideComposition.HEPG2.3.DNA.tsv.gz
    │           │   └── BCNucleotideComposition.HEPG2.3.RNA.tsv.gz
    │           ├── counts.filtered.tsv
    │           ├── counts.freqUMIs.HEPG2.1.DNA.txt
    │           ├── counts.freqUMIs.HEPG2.1.RNA.txt
    │           ├── counts.freqUMIs.HEPG2.2.DNA.txt
    │           ├── counts.freqUMIs.HEPG2.2.RNA.txt
    │           ├── counts.freqUMIs.HEPG2.3.DNA.txt
    │           ├── counts.freqUMIs.HEPG2.3.RNA.txt
    │           ├── counts.raw.tsv
    │           ├── statistic_assigned_bc_correlation_merged.fromWorkflow.default.tsv
    │           ├── statistic_assigned_counts_merged.fromWorkflow.default.tsv
    │           ├── statistic_assigned_counts_single.fromWorkflow.default.tsv
    │           ├── statistic_bc_correlation_merged.default.tsv
    │           └── statistic_oligo_correlation_merged.fromWorkflow.default.tsv
    └── logs
