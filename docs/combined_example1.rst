.. _Combined example:

.. role:: bash(code)
   :language: bash

============================
Combined workflow
============================

This example runs the assignment and the experiment/count workflow on 5'/5' WT MRPA data in the HEPG2 cell line from `Klein J., Agarwal, V., Keith, A., et al. 2019 <https://www.biorxiv.org/content/10.1101/576405v1.full.pdf>`_.

Prerequirements
======================

This example depends on the following data and software:


Installation of MPRAsnakeflow
----------------------------------------

Please install conda, the MPRAsnakeflow environment and clone the actual MPRAsnakeflow master branch. You will find more help under :ref:`Installation`.

Meta Data
___________

It is necessary to get the ordered oligo array so that each enhancer sequence can be labeled in the analysis and to trim any adaptors still in the sequence, in this case we trim off 15bp from the end of each sequence

.. code-block:: bash

    mkdir -p combined_basic/data
    cd combined_basic/data
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4237nnn/GSM4237954/suppl/GSM4237954_9MPRA_elements.fa.gz

    zcat GSM4237954_9MPRA_elements.fa.gz |awk '{ count+=1; if (count == 1) { print } else { print substr($1,1,171)}; if (count == 2) { count=0 } }' > design.fa

Reads
----------

There is one set of association sequencing for this data, which contains a forward (CRS-forward), reverse (CRS-reverse), and index (barcode) read for DNA and RNA. These data must be downloaded. All data is publically available on the short read archive (SRA). We will use SRA-toolkit to obtain the data.

.. note:: You need 10 GB disk space to download the data!

.. code-block:: bash

    conda install sra-tools
    cd combined_basic/data
    fastq-dump --gzip --split-files SRR10800986
    cd ..

For large files and unstable internet connection we reccommend the comand `prefetch` from SRA tools before running `fastq-dump`. This command is much smarter in warnings when something went wrong.

.. code-block:: bash

    conda install sra-tools
    cd combined_basic/data
    prefetch SRR10800986
    fastq-dump --gzip --split-files SRR10800986
    cd ..

.. note:: Please be sure that all files are downloaded completely without errors! Depending on your internet connection this can take a while. If you just want some data to run MPRsnakeAflow you can just limit yourself to one condition and/or just one replicate. 


With

.. code-block:: bash

    tree data


the folder should look like this:

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
    ├── SRR10800986_1.fastq.gz
    ├── SRR10800986_2.fastq.gz
    └── SRR10800986_3.fastq.gz


MPRAsnakeflow
=================================

Now we are ready to run MPRAsnakeflow and create CRS-barcode mappings and counts.

Run snakemake
------------------------------

Now we have everything at hand to run the count MPRAsnakeflow pipeline. We will run the pipeline directly in the :code:`combined_basic` folder. The MPRAsnakeflow workflow can be in a different directory. Let's assume :code:`/home/user/MPRAsnakeflow`. 

First we have to configure the config file:

.. literalinclude:: ../resources/combined_basic/config.yml
   :language: yaml


First we do a try run using snakemake :code:`-n` option. The MPRAsnakeflow command is:

.. code-block:: bash

    cd combined_basic
    conda activate mprasnakeflow
    snakemake -c 1 --use-conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile /home/user/MPRAsnakeflow/resources/combined_basic/config.yml -n

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
    assignment_getBCs                                                   1              1              1
    assignment_getInputs                                                3              1              1
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
    statistic_assigned_counts_combine_BC_assignment_stats               1              1              1
    statistic_assigned_counts_combine_BC_assignment_stats_helper        1              1              1
    statistic_assigned_counts_combine_stats_dna_rna_merge               1              1              1
    statistic_assigned_counts_combine_stats_dna_rna_merge_all           1              1              1
    statistic_bc_overlap_combine_assigned_counts                        1              1              1
    statistic_bc_overlap_combine_counts                                 1              1              1
    statistic_bc_overlap_run                                            4              1              1
    statistic_correlation_bc_counts                                     2              1              1
    statistic_correlation_calculate                                     1              1              1
    statistic_correlation_combine_bc_assigned                           1              1              1
    statistic_correlation_combine_bc_raw                                1              1              1
    statistic_correlation_combine_oligo                                 1              1              1
    statistic_counts_BC_in_RNA_DNA                                      6              1              1
    statistic_counts_BC_in_RNA_DNA_merge                                2              1              1
    statistic_counts_barcode_base_composition                           6              1              1
    statistic_counts_final                                              2              1              1
    statistic_counts_frequent_umis                                      6              1              1
    statistic_counts_stats_merge                                        2              1              1
    statistic_counts_table                                             12              1              1
    total                                                             139              1             10


When dry-drun does not give any errors we will run the workflow. We use a machine with 30 threads/cores to run the workflow. Therefore :code:`split_number` is set to 30 to parallize the workflow. Also we are using 10 threads for mapping (bwa mem). But snakemake takes care that no more than 30 threads are used.

.. code-block:: bash

    snakemake -c 30 --use-conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile /home/user/MPRAsnakeflow/resources/combined_basic/config.yml


.. note:: Please modify your code when running in a cluster environment. We have an example SLURM config file here :code:`config/sbatch.yml`.


Results
-----------------

All needed output files will be in the :code:`results/assignemnts/assocBasic` folder. The final assignment is in :code:`results/assignment/assocBasic/assignment_barcodes.exampleConfigTrueMatches.sorted.tsv.gz` or :code:`results/assignment/assocBasic/assignment_barcodes.exampleConfig.sorted.tsv.gz` depeding on the filtering in the config file. 

.. note:: Please note that for the experiment/count workflow you have to remove ambigous BCs. Therefore the file :code:`results/assignment/assocBasic/assignment_barcodes.exampleConfigTrueMatches.sorted.tsv.gz` is the correct wone


Total file tree of the results folder:

.. code-block:: text

    TODO
