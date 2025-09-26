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
    snakemake -c 1 --use-conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile /home/user/MPRAsnakeflow/resources/combined_basic/config.yml -n -q

You should see a list of rules that will be executed. This is the summary:

.. code-block:: text

    Job stats:
    job                                                             count    min threads    max threads
    ------------------------------------------------------------  -------  -------------  -------------
    all                                                                 1              1              1
    assigned_counts_assignBarcodes                                      6              1              1
    assigned_counts_dna_rna_merge                                       3              1              1
    assigned_counts_filterAssignment                                    1              1              1
    assigned_counts_make_master_tables                                  1              1              1
    assignment_bwa_ref                                                  1              1              1
    assignment_fastq_split                                              3              1              1
    assignment_filter                                                   1              1              1
    assignment_flagstat                                                 1              1              1
    assignment_mapping_bwa_getBCs                                       1              1              1
    assignment_idx_bam                                                  1              1              1
    assignment_mapping                                                  1             10             10
    assignment_merge                                                   30              1              1
    assignment_statistic_assignedCounts                                 1              1              1
    assignment_statistic_assignment                                     1              1              1
    assignment_statistic_totalCounts                                    1              1              1
    counts_create_BAM_umi                                               6              1              1
    counts_dna_rna_merge_counts                                         6              1              1
    counts_filter_counts                                                6              1              1
    counts_final_counts_umi                                             6              1              1
    counts_raw_counts_umi                                               6              1              1
    statistic_counts_frequent_umis                                      6              1              1
    total                                                             136              1             10

When the dry-run does not give any errors, run the workflow. We use a machine with 30 threads/cores to run the workflow. Therefore, :code:`split_number` is set to 30 to parallelize the workflow. Also, we are using 10 threads for mapping (bwa mem). Snakemake ensures that no more than 30 threads are used.

.. code-block:: bash

    snakemake -c 30 --use-conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile /home/user/MPRAsnakeflow/resources/combined_basic/config.yml

.. note:: Please modify your code when running in a cluster environment. We have an example SLURM config file here: :code:`config/sbatch.yml`.

Results
-------
All needed output files will be in the :code:`results/assignment/assocBasic` folder for assignment results. The folder :code:`results/experiments/countBasic` contains the count results.

To generate a final report, use the following command:

.. code-block:: bash

    snakemake --config config.yml --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --report report.html

This HTML report contains information about the Snakemake run and integrates statistics tables and plots.

Total file tree of the results folder:

.. code-block:: text

    results
    ├── assignment
    │   └── assocBasic
    │       ├── aligned_merged_reads.bam
    │       ├── aligned_merged_reads.bam.bai
    │       ├── assignment_barcodes.default.tsv.gz
    │       ├── barcodes_incl_other.tsv.gz
    │       ├── reference
    │       │   ├── reference.fa
    │       │   ├── reference.fa.amb
    │       │   ├── reference.fa.ann
    │       │   ├── reference.fa.bwt
    │       │   ├── reference.fa.dict
    │       │   ├── reference.fa.fai
    │       │   ├── reference.fa.pac
    │       │   └── reference.fa.sa
    │       └── statistic
    │           ├── assigned_counts.default.tsv
    │           ├── assignment
    │           │   └── bam_stats.txt
    │           ├── assignment.default.png
    │           ├── assignment.default.tsv.gz
    │           └── total_counts.tsv.gz
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
