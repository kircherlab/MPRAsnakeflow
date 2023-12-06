.. _Count example:

.. role:: bash(code)
      :language: bash

=========================
Basic Experiment workflow
=========================

This example runs the count workflow on 5'/5' WT MPRA data in the HepG2 cell line from `Klein J., Agarwal, V., Keith, A., et al. 2019 <https://www.biorxiv.org/content/10.1101/576405v1.full.pdf>`_.

Prerequirements
======================

This example depends on the following data and software:


Installing MPRAsnakeflow
------------------------

Please install conda, the MPRAsnakeflow environment, and clone the actual ``MPRAsnakeflow`` master branch. You will find more help under :ref:`Installation`.

Producing an association (.tsv.gz) file 
----------------------------------------
This workflow requires a python dictionary of candidate regulatory sequence (CRS) mapped to their barcodes in a tab separated (.tsv) format. For this example the file can be generated using :ref:`Assignment example` or it can be found in :code:`resources/count_basic` folder in `MPRAsnakelfow <https://github.com/kircherlab/MPRAsnakeflow/>`_(file :code:`SRR10800986_barcodes_to_coords.tsv.gz`).

Alternatively, if the association file is in pickle (.pickle) format because you used MPRAflow, you can convert the same file to .tsv.gz format with the in-built function in MPRsnakeflow with the following code:

.. code-block:: bash
    
    conda activate mprasnakeflow
    python assignment_pickle_to_tsv.py --input <assignment_file>.pickle | sort | uniq | gzip -c > <assignment_file>.tsv.gz


Design (.fa) file
-----------------

    File can be generated using the :ref:`Assignment example` or downloaded from the :code:`resources/count_basic` folder in `MPRAsnakelfow <https://github.com/kircherlab/MPRAsnakeflow/>`_.



Reads
-----

There is one condition (HEPG2) with three technical replicates. Each replicate contains a forward (barcode-forward), reverse (barcode-reverse), and index (unique molecular identifier) read for DNA and RNA. These data must be downloaded. All data is publically available on the short read archive (SRA). We will use SRA-toolkit to obtain the data.

.. note:: You need 9 GB disk space to download the data and upwards of 50 GB to proccess it!

.. code-block:: bash

    conda install sra-tools
    mkdir -p count_basic/data
    cd count_basic/data
    fastq-dump --gzip --split-files SRR10800881 SRR10800882 SRR10800883 SRR10800884 SRR10800885 SRR10800886
    cd ..

For large files and unstable internet connection we reccommend the comand ``prefetch`` from SRA tools before running ``fastq-dump``. This command is much smarter in warnings when something went wrong.

.. code-block:: bash

    conda install sra-tools
    cd count_basic/data
    prefetch SRR10800881 SRR10800882 SRR10800883 SRR10800884 SRR10800885 SRR10800886
    fastq-dump --gzip --split-files SRR10800986
    cd ..


.. note:: Please be sure that all files are downloaded completely without errors! Depending on your internet connection this can take a while. If you just want some data to run MPRAsnakeflow, you can just limit yourself to one condition and/or just one replicate.

The data folder view can be seen with the following command:

.. code-block:: bash

    tree data


The folder should look like this:

.. code-block:: text

    data
    ├── design.fa
    ├── SRR10800881_1.fastq.gz
    ├── SRR10800881_2.fastq.gz
    ├── SRR10800881_3.fastq.gz
    ├── SRR10800882_1.fastq.gz
    ├── SRR10800882_2.fastq.gz
    ├── SRR10800882_3.fastq.gz
    ├── SRR10800883_1.fastq.gz
    ├── SRR10800883_2.fastq.gz
    ├── SRR10800883_3.fastq.gz
    ├── SRR10800884_1.fastq.gz
    ├── SRR10800884_2.fastq.gz
    ├── SRR10800884_3.fastq.gz
    ├── SRR10800885_1.fastq.gz
    ├── SRR10800885_2.fastq.gz
    ├── SRR10800885_3.fastq.gz
    ├── SRR10800886_1.fastq.gz
    ├── SRR10800886_2.fastq.gz
    ├── SRR10800886_3.fastq.gz
    └── SRR10800986_barcodes_to_coords.tsv.gz

Here is an overview of the files:

.. csv-table:: HEPG2 data
   :header: "Condition", "GEO Accession", "SRA Accession", SRA Runs
   :widths: 40, 10, 10, 20

   "HEPG2-DNA-1: HEPG2 DNA replicate 1", GSM4237863, SRX7474781, "SRR10800881"
   "HEPG2-RNA-1: HEPG2 RNA replicate 1", GSM4237864, SRX7474782, "SRR10800882"
   "HEPG2-DNA-2: HEPG2 DNA replicate 2", GSM4237865, SRX7474783, "SRR10800883"
   "HEPG2-RNA-2: HEPG2 RNA replicate 2", GSM4237866, SRX7474784, "SRR10800884"
   "HEPG2-DNA-3: HEPG2 DNA replicate 3", GSM4237867, SRX7474785, "SRR10800885"
   "HEPG2-RNA-3: HEPG2 RNA replicate 3", GSM4237868, SRX7474786, "SRR10800886"



Run MPRAsnakeflow
=================

Now we are close to starting MPRAsnakeflow and count the number of barcodes. But before we need to generate an environment (.csv) file to tell snakemake the conditions, replicates and the corresponding reads.

Creating experiment.csv
---------------------------

Our experiment file looks exactly like this:

.. code-block:: text

    Condition,Replicate,DNA_BC_F,DNA_UMI,DNA_BC_R,RNA_BC_F,RNA_UMI,RNA_BC_R
    HEPG2,1,SRR10800881_1.fastq.gz,SRR10800881_2.fastq.gz,SRR10800881_3.fastq.gz,SRR10800882_1.fastq.gz,SRR10800882_2.fastq.gz,SRR10800882_3.fastq.gz
    HEPG2,2,SRR10800883_1.fastq.gz,SRR10800883_2.fastq.gz,SRR10800883_3.fastq.gz,SRR10800884_1.fastq.gz,SRR10800884_2.fastq.gz,SRR10800884_3.fastq.gz
    HEPG2,3,SRR10800885_1.fastq.gz,SRR10800885_2.fastq.gz,SRR10800885_3.fastq.gz,SRR10800886_1.fastq.gz,SRR10800886_2.fastq.gz,SRR10800886_3.fastq.gz

Save it into the :code:`count_basic/data` folder under :code:`experiment.csv`.

MPRAsnakeflow
=================================

Now we have everything at hand to run the count MPRAsnakeflow pipeline. We will run the pipeline directly in the :code:`count_basic` folder. The MPRAsnakeflow workflow can be in a different directory. Let's assume :code:`/home/user/MPRAsnakeflow`.  The MPRAsnakeflow count command is:


First we have to configure the config file and save it to the :code:`count_basic` folder. The config file is a simple text file with the following content:

.. literalinclude:: ../resources/count_basic/config.yml
   :language: yaml


First we do a try run using snakemake :code:`-n` option. The MPRAsnakeflow command is:

.. code-block:: bash

    cd count_basic
    conda activate mprasnakeflow
    snakemake -c 1 --use-conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile /home/user/MPRAsnakeflow/resources/count_basic/config.yml -n -q

You should see a list of rules that will be executed. This is the summary:

.. code-block:: text

    Building DAG of jobs...
    Job stats:
    job                                                             count
    ------------------------------------------------------------  -------
    all                                                                 1
    assigned_counts_assignBarcodes                                      6
    assigned_counts_combine_replicates                                  2
    assigned_counts_combine_replicates_barcode_output                   1
    assigned_counts_dna_rna_merge                                       3
    assigned_counts_filterAssignment                                    1
    assigned_counts_make_master_tables                                  1
    counts_dna_rna_merge_counts                                         6
    counts_filter_counts                                                6
    counts_final_counts                                                 6
    counts_umi_create_BAM                                               6
    counts_umi_raw_counts                                               6
    statistic_assigned_counts_combine_BC_assignment_stats               1
    statistic_assigned_counts_combine_BC_assignment_stats_helper        1
    statistic_assigned_counts_combine_stats_dna_rna_merge               1
    statistic_assigned_counts_combine_stats_dna_rna_merge_all           1
    statistic_bc_overlap_combine_assigned_counts                        1
    statistic_bc_overlap_combine_counts                                 1
    statistic_bc_overlap_run                                            4
    statistic_correlation_bc_counts                                     2
    statistic_correlation_bc_counts_hist                                2
    statistic_correlation_calculate                                     1
    statistic_correlation_combine_bc_assigned                           1
    statistic_correlation_combine_bc_raw                                1
    statistic_correlation_combine_oligo                                 1
    statistic_correlation_hist_box_plots                                1
    statistic_counts_BC_in_RNA_DNA                                      6
    statistic_counts_BC_in_RNA_DNA_merge                                2
    statistic_counts_barcode_base_composition                           6
    statistic_counts_final                                              2
    statistic_counts_frequent_umis                                      6
    statistic_counts_stats_merge                                        2
    statistic_counts_table                                             12
    total                                                             100

When dry-drun does not give any errors we will run the workflow. We use a machine with 30 threads/cores to run the workflow. The MPRAsnakeflow command is:

.. code-block:: bash

    snakemake -c 30 --use-conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile /home/user/MPRAsnakeflow/resources/count_basic/config.yml

.. note:: Please modify your code when running in a cluster environment. We have an example SLURM config file here :code:`config/sbatch.yml`.

If everything works fine the 29 rules showed above will run. Everything starting with :code:`counts_` beolngs to raw count rules, with :code:`assigned_counts_` to counts assigned to the assignment and :code:`statistic_` to statistics. Here is a brief description of the rules.

all
    The overall all rule. Here is defined what final output files are expected.
counts_create_BAM_umi
    Create a BAM file from FASTQ input, merge FW and REV read and save UMI in XI flag.
counts_raw_counts_umi
    Counting BCsxUMIs from the BAM files.
counts_filter_counts
    Filter the counts to BCs only of the correct length (defined in the config file).
counts_final_counts_umi
    Discarding PCR duplicates (taking BCxUMI only one time). Final result of counts can be found here: :code:`results/experiments/exampleCount/counts/HepG2_<1,2,3>_<DNA/RNA>_filtered_counts.tsv.gz`.
counts_dna_rna_merge_counts
    Merge DNA and RNA counts together.
    This is done in two ways. First no not allow zeros in DNA or RNA BCs (when :code:`min_counts` is not zero for DNA and RNA).
    Second with zeros, so a BC can be defined only in the DNA or RNA (when :code:`min_counts` is zero for DNA or RNA)
assigned_counts_filterAssignment
    Use only unique assignments.
assigned_counts_combine_replicates
    TODO
assigned_counts_combine_replicates_barcode_output
    TODO
assigned_counts_assignBarcodes
    Assign RNA and DNA barcodes seperately to make the statistic for assigned.
assigned_counts_dna_rna_merge
    Assign merged RNA/DNA barcodes. Filter BC depending on the min_counts option. Output for each replicate is here: :code:`results/experiments/exampleCount/assigned_counts/fromFile/default/HepG2_<1,2,3>_merged_assigned_counts.tsv.gz`.
assigned_counts_make_master_tables
    Final master table with all replicates combined. Output is here: :code:`results/experiments/exampleCount/assigned_counts/fromFile/default/HepG2_allreps_merged.tsv.gz` and using the :code:`bc-threshold` here :code:`results/experiments/exampleCount/assigned_counts/fromFile/default/HepG2_allreps_minThreshold_merged.tsv.gz`.
statistic_assigned_counts_combine_BC_assignment_stats
    TODO
statistic_assigned_counts_combine_BC_assignment_stats_helper
    TODO
statistic_assigned_counts_combine_stats_dna_rna_merge
    TODO
statistic_assigned_counts_combine_stats_dna_rna_merge_all
    TODO
statistic_bc_overlap_combine_assigned_counts
    TODO
statistic_bc_overlap_combine_counts
    TODO
statistic_bc_overlap_run
    TODO
statistic_correlation_bc_counts
    TODO
statistic_correlation_calculate
    TODO
statistic_correlation_combine_bc_assigned
    TODO
statistic_correlation_combine_bc_raw
    TODO
statistic_correlation_combine_oligo
    TODO
statistic_counts_BC_in_RNA_DNA
    TODO
statistic_counts_BC_in_RNA_DNA_merge
    TODO
statistic_counts_barcode_base_composition
    TODO
statistic_counts_final
    TODO
statistic_counts_frequent_umis
    TODO
statistic_counts_stats_merge
    TODO
statistic_counts_table
    TODO




.. todo:: Rules not correct in example experiment workflow

Results
-----------------

All needed output files will be in the :code:`results/assignment/countBasic`folder.


To generate a final report, the following code can be used

.. code-block:: bash

    snakemake --config config.yml --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --report report.html 

This html contains als information about the snakemake run and integrates statistic tables and plots.

Total file tree of the results folder:

.. code-block:: text
    results
    └── experiments
        └── countBasic
            ├── assigned_counts
            │   └── fromWorkflow
            │       ├── default
            │       │   ├── HEPG2_1_merged_assigned_counts.tsv.gz
            │       │   ├── HEPG2_2_merged_assigned_counts.tsv.gz
            │       │   ├── HEPG2_3_merged_assigned_counts.tsv.gz
            │       │   ├── HEPG2_allreps_merged.tsv.gz
            │       │   └── HEPG2_allreps_minThreshold_merged.tsv.gz
            │       ├── HEPG2_1_DNA_final_counts.config.default.tsv.gz
            │       ├── HEPG2_1.merged.config.default.tsv.gz
            │       ├── HEPG2_1_RNA_final_counts.config.default.tsv.gz
            │       ├── HEPG2_2_DNA_final_counts.config.default.tsv.gz
            │       ├── HEPG2_2.merged.config.default.tsv.gz
            │       ├── HEPG2_2_RNA_final_counts.config.default.tsv.gz
            │       ├── HEPG2_3_DNA_final_counts.config.default.tsv.gz
            │       ├── HEPG2_3.merged.config.default.tsv.gz
            │       └── HEPG2_3_RNA_final_counts.config.default.tsv.gz
            ├── assignment
            │   └── fromWorkflow.tsv.gz
            ├── counts
            │   ├── HEPG2_1_DNA.bam
            │   ├── HEPG2_1_DNA_filtered_counts.tsv.gz
            │   ├── HEPG2_1_DNA_final_counts.tsv.gz
            │   ├── HEPG2_1_DNA_raw_counts.tsv.gz
            │   ├── HEPG2_1.merged.config.default.tsv.gz
            │   ├── HEPG2_1_RNA.bam
            │   ├── HEPG2_1_RNA_filtered_counts.tsv.gz
            │   ├── HEPG2_1_RNA_final_counts.tsv.gz
            │   ├── HEPG2_1_RNA_raw_counts.tsv.gz
            │   ├── HEPG2_2_DNA.bam
            │   ├── HEPG2_2_DNA_filtered_counts.tsv.gz
            │   ├── HEPG2_2_DNA_final_counts.tsv.gz
            │   ├── HEPG2_2_DNA_raw_counts.tsv.gz
            │   ├── HEPG2_2.merged.config.default.tsv.gz
            │   ├── HEPG2_2_RNA.bam
            │   ├── HEPG2_2_RNA_filtered_counts.tsv.gz
            │   ├── HEPG2_2_RNA_final_counts.tsv.gz
            │   ├── HEPG2_2_RNA_raw_counts.tsv.gz
            │   ├── HEPG2_3_DNA.bam
            │   ├── HEPG2_3_DNA_filtered_counts.tsv.gz
            │   ├── HEPG2_3_DNA_final_counts.tsv.gz
            │   ├── HEPG2_3_DNA_raw_counts.tsv.gz
            │   ├── HEPG2_3.merged.config.default.tsv.gz
            │   ├── HEPG2_3_RNA.bam
            │   ├── HEPG2_3_RNA_filtered_counts.tsv.gz
            │   ├── HEPG2_3_RNA_final_counts.tsv.gz
            │   └── HEPG2_3_RNA_raw_counts.tsv.gz
            └── statistic
                ├── assigned_counts
                │   └── fromWorkflow
                │       ├── default
                │       │   ├── combined
                │       │   │   └── HEPG2_merged_assigned_counts.statistic.tsv.gz
                │       │   ├── HEPG2_1_merged_assigned_counts.statistic.tsv.gz
                │       │   ├── HEPG2_2_merged_assigned_counts.statistic.tsv.gz
                │       │   ├── HEPG2_3_merged_assigned_counts.statistic.tsv.gz
                │       │   ├── HEPG2_all_barcodesPerInsert_box_minThreshold.png
                │       │   ├── HEPG2_all_barcodesPerInsert_box.png
                │       │   ├── HEPG2_average_allreps_merged.tsv.gz
                │       │   ├── HEPG2_barcodesPerInsert.png
                │       │   ├── HEPG2_correlation_minThreshold.tsv
                │       │   ├── HEPG2_correlation.tsv
                │       │   ├── HEPG2_DNA_pairwise_minThreshold.png
                │       │   ├── HEPG2_DNA_pairwise.png
                │       │   ├── HEPG2_group_barcodesPerInsert_box_minThreshold.png
                │       │   ├── HEPG2_group_barcodesPerInsert_box.png
                │       │   ├── HEPG2_Ratio_pairwise_minThreshold.png
                │       │   ├── HEPG2_Ratio_pairwise.png
                │       │   ├── HEPG2_RNA_pairwise_minThreshold.png
                │       │   └── HEPG2_RNA_pairwise.png
                │       ├── HEPG2_1_DNA_default.statistic.tsv.gz
                │       ├── HEPG2_1_RNA_default.statistic.tsv.gz
                │       ├── HEPG2_2_DNA_default.statistic.tsv.gz
                │       ├── HEPG2_2_RNA_default.statistic.tsv.gz
                │       ├── HEPG2_3_DNA_default.statistic.tsv.gz
                │       └── HEPG2_3_RNA_default.statistic.tsv.gz
                ├── barcode
                │   ├── assigned_counts
                │   │   └── fromWorkflow
                │   │       ├── HEPG2_default_barcode_correlation.tsv
                │   │       ├── HEPG2_default_barcode_DNA_pairwise.png
                │   │       ├── HEPG2_default_barcode_Ratio_pairwise.png
                │   │       ├── HEPG2_default_barcode_RNA_pairwise.png
                │   │       ├── HEPG2_default_DNA_perBarcode.png
                │   │       └── HEPG2_default_RNA_perBarcode.png
                │   └── counts
                │       ├── HEPG2_default_barcode_correlation.tsv
                │       ├── HEPG2_default_barcode_DNA_pairwise.png
                │       ├── HEPG2_default_barcode_Ratio_pairwise.png
                │       ├── HEPG2_default_barcode_RNA_pairwise.png
                │       ├── HEPG2_default_DNA_perBarcode.png
                │       └── HEPG2_default_RNA_perBarcode.png
                ├── bc_overlap.assigned_counts.default.fromWorkflow.tsv
                ├── bc_overlap.counts.default.tsv
                ├── counts
                │   ├── BCNucleotideComposition.HEPG2_1_DNA.tsv.gz
                │   ├── BCNucleotideComposition.HEPG2_1_RNA.tsv.gz
                │   ├── BCNucleotideComposition.HEPG2_2_DNA.tsv.gz
                │   ├── BCNucleotideComposition.HEPG2_2_RNA.tsv.gz
                │   ├── BCNucleotideComposition.HEPG2_3_DNA.tsv.gz
                │   └── BCNucleotideComposition.HEPG2_3_RNA.tsv.gz
                ├── counts.filtered.tsv
                ├── counts.freqUMIs.HEPG2_1_DNA.txt
                ├── counts.freqUMIs.HEPG2_1_RNA.txt
                ├── counts.freqUMIs.HEPG2_2_DNA.txt
                ├── counts.freqUMIs.HEPG2_2_RNA.txt
                ├── counts.freqUMIs.HEPG2_3_DNA.txt
                ├── counts.freqUMIs.HEPG2_3_RNA.txt
                ├── counts.raw.tsv
                ├── statistic_assigned_bc_correlation_merged_fromWorkflow_default.tsv
                ├── statistic_assigned_counts_merged_fromWorkflow_default.tsv
                ├── statistic_assigned_counts_single_fromWorkflow_default.tsv
                ├── statistic_bc_correlation_merged_default.tsv
                └── statistic_oligo_correlation_merged_fromWorkflow_default.tsv
