.. _Count example:

.. role:: bash(code)
      :language: bash

=========================
Basic Experiment Workflow
=========================

This example runs the count workflow on 5'/5' WT MPRA data in the HepG2 cell line from `Klein J., Agarwal, V., Keith, A., et al. 2019 <https://www.biorxiv.org/content/10.1101/576405v1.full.pdf>`_.

Prerequisites
=============

This example depends on the following data and software:

Installing MPRAsnakeflow
------------------------
Please install conda, the MPRAsnakeflow environment, and clone the actual ``MPRAsnakeflow`` master branch. You will find more help under :ref:`Installation`.

Producing an Association (.tsv.gz) File
---------------------------------------
This workflow requires a Python dictionary of candidate regulatory sequences (CRS) mapped to their barcodes in a tab-separated (.tsv.gz) format. For this example, the file can be generated using :ref:`Assignment example` or found in the :code:`resources/count_basic` folder in `MPRAsnakeflow <https://github.com/kircherlab/MPRAsnakeflow/>`_ (file: :code:`SRR10800986_barcodes_to_coords.tsv.gz`).

Alternatively, if the association file is in pickle (.pickle) format because you used MPRAflow, you can convert it to .tsv.gz format using the in-built function in MPRAsnakeflow with the following code:

.. code-block:: bash

    conda activate mprasnakeflow
    python assignment_pickle_to_tsv.py --input <assignment_file>.pickle | sort | uniq | gzip -c > <assignment_file>.tsv.gz

Design (.fa) File
-----------------
The design file can be generated using the :ref:`Assignment example` or downloaded from the :code:`resources/count_basic` folder in `MPRAsnakeflow <https://github.com/kircherlab/MPRAsnakeflow/>`_.

Reads
-----
There is one condition (HepG2) with three technical replicates. Each replicate contains a forward (barcode-forward), reverse (barcode-reverse), and index (unique molecular identifier) read for DNA and RNA. These data must be downloaded. All data is publicly available on the Short Read Archive (SRA). We will use SRA-toolkit to obtain the data.

.. note:: You need 9 GB of disk space to download the data and upwards of 50 GB to process it!

.. code-block:: bash

    conda install sra-tools
    mkdir -p count_basic/data
    cd count_basic/data
    fastq-dump --gzip --split-files SRR10800881 SRR10800882 SRR10800883 SRR10800884 SRR10800885 SRR10800886
    cd ..

For large files and unstable internet connections, we recommend the command ``prefetch`` from SRA tools before running ``fastq-dump``. This command provides better warnings when something goes wrong.

.. code-block:: bash

    conda install sra-tools
    cd count_basic/data
    prefetch SRR10800881 SRR10800882 SRR10800883 SRR10800884 SRR10800885 SRR10800886
    fastq-dump --gzip --split-files SRR10800986
    cd ..

.. note:: Please ensure that all files are downloaded completely without errors! Depending on your internet connection, this can take a while. If you just want some data to run MPRAsnakeflow, you can limit yourself to one condition and/or just one replicate.

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

.. csv-table:: HepG2 Data
   :header: "Condition", "GEO Accession", "SRA Accession", "SRA Runs"
   :widths: 40, 10, 10, 20

   "HepG2-DNA-1: HepG2 DNA replicate 1", GSM4237863, SRX7474781, "SRR10800881"
   "HepG2-RNA-1: HepG2 RNA replicate 1", GSM4237864, SRX7474782, "SRR10800882"
   "HepG2-DNA-2: HepG2 DNA replicate 2", GSM4237865, SRX7474783, "SRR10800883"
   "HepG2-RNA-2: HepG2 RNA replicate 2", GSM4237866, SRX7474784, "SRR10800884"
   "HepG2-DNA-3: HepG2 DNA replicate 3", GSM4237867, SRX7474785, "SRR10800885"
   "HepG2-RNA-3: HepG2 RNA replicate 3", GSM4237868, SRX7474786, "SRR10800886"

Run MPRAsnakeflow
=================
Now we are ready to start MPRAsnakeflow and count the number of barcodes. First, we need to generate an experiment file (`experiment.csv`) to specify the conditions, replicates, and corresponding reads.

Creating experiment.csv
-----------------------
Our experiment file looks like this:

.. code-block:: text

    Condition,Replicate,DNA_BC_F,DNA_UMI,DNA_BC_R,RNA_BC_F,RNA_UMI,RNA_BC_R
    HepG2,1,SRR10800881_1.fastq.gz,SRR10800881_2.fastq.gz,SRR10800881_3.fastq.gz,SRR10800882_1.fastq.gz,SRR10800882_2.fastq.gz,SRR10800882_3.fastq.gz
    HepG2,2,SRR10800883_1.fastq.gz,SRR10800883_2.fastq.gz,SRR10800883_3.fastq.gz,SRR10800884_1.fastq.gz,SRR10800884_2.fastq.gz,SRR10800884_3.fastq.gz
    HepG2,3,SRR10800885_1.fastq.gz,SRR10800885_2.fastq.gz,SRR10800885_3.fastq.gz,SRR10800886_1.fastq.gz,SRR10800886_2.fastq.gz,SRR10800886_3.fastq.gz

Save it into the :code:`count_basic/data` folder as :code:`experiment.csv`.

MPRAsnakeflow
=============
Now we have everything ready to run the count workflow of MPRAsnakeflow. Run the pipeline directly in the :code:`count_basic` folder. The MPRAsnakeflow workflow can be in a different directory. Let's assume :code:`/home/user/MPRAsnakeflow`.

First, configure the config file and save it to the :code:`count_basic` folder. The config file is a simple text file with the following content:

.. literalinclude:: ../resources/count_basic/config.yml
   :language: yaml

Perform a dry-run using Snakemake's :code:`-n` option:

.. code-block:: bash

    cd count_basic
    conda activate mprasnakeflow
    snakemake -c 1 --use-conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile /home/user/MPRAsnakeflow/resources/count_basic/config.yml -n -q

If the dry-run does not give any errors, run the workflow using 30 threads:

.. code-block:: bash

    snakemake -c 30 --use-conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile /home/user/MPRAsnakeflow/resources/count_basic/config.yml

.. note:: Please modify your code when running in a cluster environment. We have an example SLURM config file here: :code:`config/sbatch.yml`.

Results
-------
All output files will be in the :code:`results/experiments/countBasic` folder.

To generate a final report, use the following command:

.. code-block:: bash

    snakemake --config config.yml --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --report report.html

This HTML report contains information about the Snakemake run and integrates statistics tables and plots.

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
