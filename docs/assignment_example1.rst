.. _Assignment example:

.. role:: bash(code)
   :language: bash

============================
Basic assignment workflow
============================

This example runs the assignment workflow on 5'/5' WT MRPA data in the HEPG2 cell line from `Klein J., Agarwal, V., Keith, A., et al. 2019 <https://www.biorxiv.org/content/10.1101/576405v1.full.pdf>`_.

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

    mkdir -p Assoc_Basic/data
    cd Assoc_Basic/data
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4237nnn/GSM4237954/suppl/GSM4237954_9MPRA_elements.fa.gz

    zcat GSM4237954_9MPRA_elements.fa.gz |awk '{ count+=1; if (count == 1) { print } else { print substr($1,1,171)}; if (count == 2) { count=0 } }' > design.fa

Reads
----------

There is one set of association sequencing for this data, which contains a forward (CRS-forward), reverse (CRS-reverse), and index (barcode) read for DNA and RNA. These data must be downloaded. All data is publically available on the short read archive (SRA). We will use SRA-toolkit to obtain the data.

.. note:: You need 10 GB disk space to download the data!

.. code-block:: bash

    conda install sra-tools
    cd Assoc_Basic/data
    fastq-dump --gzip --split-files SRR10800986
    cd ..

For large files and unstable internet connection we reccommend the comand `prefetch` from SRA tools before running `fastq-dump`. This command is much smarter in warnings when something went wrong.

.. code-block:: bash

    conda install sra-tools
    cd Assoc_Basic/data
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

Here is an overview of the files:

.. csv-table:: HEPG2 association data
   :header: "Condition", "GEO Accession", "SRA Accession", SRA Runs
   :widths: 40, 10, 10, 20

   "HEPG2-association: HEPG2 library association", GSM4237954, SRX7474872, "SRR10800986"


MPRAsnakeflow
=================================

Now we are ready to run MPRAsnakeflow and create CRS-barcode mappings.

Run snakemake
------------------------------

Now we have everything at hand to run the count MPRAsnakeflow pipeline. We will run the pipeline directly in the :code:`Assoc_Basic` folder. The MPRAsnakeflow workflow can be in a different directory. Let's assume :code:`/home/user/MPRAsnakeflow`. 

First we have to configure the config file:

.. literalinclude:: ../resources/assoc_basic/config.yml
   :language: yaml


First we do a try run using snakemake :code:`-n` option. The MPRAsnakeflow command is:

.. code-block:: bash

    cd Assoc_Basic
    conda activate mprasnakeflow
    snakemake -c 1 --use-conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile /home/user/MPRAsnakeflow/resources/assoc_basic/config.yml -n

You should see a list of rules that will be executed. This is the summary:

.. code-block:: text
   
   Job stats:
   job                                    count    min threads    max threads
   -----------------------------------  -------  -------------  -------------
   all                                        1              1              1
   assignment_bwa_ref                         1              1              1
   assignment_fastq_split                     3              1              1
   assignment_filter                          2              1              1
   assignment_flagstat                        1              1              1
   assignment_getBCs                          1              1              1
   assignment_getInputs                       3              1              1
   assignment_idx_bam                         1              1              1
   assignment_mapping                         1              1              1
   assignment_merge                           30             10             10
   assignment_statistic_assignedCounts        2              1              1
   assignment_statistic_assignment            2              1              1
   assignment_statistic_totalCounts           1              1              1
   total                                     49              1              1


When dry-drun does not give any errors we will run the workflow. We use a machine with 30 threads/cores to run the workflow. Therefore :code:`split_number` is set to 30 to parallize the workflow. Also we are using 10 threads for mapping (bwa mem). But snakemake takes care that no more than 30 threads are used.

.. code-block:: bash

    snakemake -c 30 --use-conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile /home/user/MPRAsnakeflow/resources/assoc_basic/config.yml


.. note:: Please modify your code when running in a cluster environment. We have an example SLURM config file here :code:`config/sbatch.yml`.

If everything works fine the 13 rules showed above will run:

all
   The overall all rule. Here is defined what final output files are expected.
assignment_bwa_ref
   Create mapping reference for BWA from design file.
assignment_fastq_split
   Split the fastq files into n files for parallelisation. N is given by split_read in the configuration file.
assignment_getInputs
   Concat the input fastq files per R1,R2,R3. If only single fastq file is provided a symbolic link is created.
assignment_merge
   Merge the FW,REV and BC fastq files into one. Extract the index sequence from the middle and end of an Illumina run. Separates reads for Paired End runs. Merge/Adapter trim reads stored in BAM.
assignment_mapping
   Map the reads to the reference.
assignment_idx_bam
   Index the BAM file
assignment_flagstat
   Run samtools flagstat. Results are in :code:`results/assignment/assoc_basic/statistic/assignment/bam_stats.txt`
assignment_getBCs
   Get the barcodes (not filtered). Results are in :code:`results/assignment/assoc_basic/barcodes_incl_other.sorted.tsv.gz`
assignment_statistic_totalCounts
   Statistic of the total (unfiltered counts). Results are in :code:`results/assignment/assoc_basic/statistic/total_counts.tsv.gz`
assignment_filter
   Filter the barcodes file based on the config given in the config-file. Results for this run are here :code:`results/assignment/assoc_basic/assignment_barcodes.example_config_true_matches.sorted.tsv.gz` (example_config_true_matches) and here :code:`results/assignment/assoc_basic/assignment_barcodes.example_config.sorted.tsv.gz` (example_config)
assignment_statistic_assignedCounts
   Statistic of filtered the assigned counts. Result is here :code:`results/assignment/assoc_basic/statistic/assigned_counts.example_config_true_matches.tsv.gz` (example_config_true_matches) or :code:`results/assignment/assoc_basic/statistic/assigned_counts.example_config.tsv.gz` (example_config)
assignment_statistic_assignment
   Statistic of the filtered assignment.  Result is here :code:`results/assignment/assoc_basic/statistic/assignment.example_config_true_matches.tsv.gz` and a plot here :code:`results/assignment/assoc_basic/statistic/assignment.example_config_true_matches.png`. (also files are available for the config :code:`example_config`).

Results
-----------------

All needed output files will be in the :code:`results/assignemnts/assoc_basic` folder. The final assignment is in :code:`results/assignment/assoc_basic/assignment_barcodes.example_config_true_matches.sorted.tsv.gz` or :code:`results/assignment/assoc_basic/assignment_barcodes.example_config.sorted.tsv.gz` depeding on the filtering in the config file. 

.. note:: Please note that for the experiment/count workflow you have to remove ambigous BCs. Therefore the file :code:`results/assignment/assoc_basic/assignment_barcodes.example_config_true_matches.sorted.tsv.gz` is the correct wone


Total file tree of the results folder:

.. code-block:: text

    ├── assignment
    │   └── assoc_basic
    │       ├── aligned_merged_reads.bam
    │       ├── aligned_merged_reads.bam.bai
    │       ├── assignment_barcodes.example_config.sorted.tsv.gz
    │       ├── assignment_barcodes.example_config_true_matches.sorted.tsv.gz
    │       ├── barcodes_incl_other.sorted.tsv.gz
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
    │           ├── assigned_counts.example_config_true_matches.tsv.gz
    │           ├── assigned_counts.example_config.tsv.gz
    │           ├── assignment
    │           │   └── bam_stats.txt
    │           ├── assignment.example_config.png
    │           ├── assignment.example_config_true_matches.png
    │           ├── assignment.example_config_true_matches.tsv.gz
    │           ├── assignment.example_config.tsv.gz
    │           └── total_counts.tsv.gz
