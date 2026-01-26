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

.. literalinclude:: ../../resources/count_basic/config.yml
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

All output files will be in the :code:`results/experiments/countBasic` folder. A nice overview (QC report) is shown in ::code:`results/experiments/countBasic/qc_report.default.html`. This HTML report contains information about statistics tables and plots. You can find an example qc report here: `Example Experiment QC report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/count_example1.qc_report.default.html>`_.

Total file tree of the results folder:

.. code-block:: text

    results/
    ├── experiments
    │   └── exampleCount
    │       ├── assigned_counts
    │       │   └── fromFile
    │       │       ├── HEPG2.1.DNA.final_counts.config.default.tsv.gz
    │       │       ├── HEPG2.1.DNA.final_counts.config.outlierZscore.tsv.gz
    │       │       ├── HEPG2.1.RNA.final_counts.config.default.tsv.gz
    │       │       ├── HEPG2.1.RNA.final_counts.config.outlierZscore.tsv.gz
    │       │       ├── HEPG2.1.merged.config.default.tsv.gz
    │       │       ├── HEPG2.1.merged.config.outlierZscore.tsv.gz
    │       │       ├── HEPG2.2.DNA.final_counts.config.default.tsv.gz
    │       │       ├── HEPG2.2.DNA.final_counts.config.outlierZscore.tsv.gz
    │       │       ├── HEPG2.2.RNA.final_counts.config.default.tsv.gz
    │       │       ├── HEPG2.2.RNA.final_counts.config.outlierZscore.tsv.gz
    │       │       ├── HEPG2.2.merged.config.default.tsv.gz
    │       │       ├── HEPG2.2.merged.config.outlierZscore.tsv.gz
    │       │       ├── HEPG2.3.DNA.final_counts.config.default.tsv.gz
    │       │       ├── HEPG2.3.DNA.final_counts.config.outlierZscore.tsv.gz
    │       │       ├── HEPG2.3.RNA.final_counts.config.default.tsv.gz
    │       │       ├── HEPG2.3.RNA.final_counts.config.outlierZscore.tsv.gz
    │       │       ├── HEPG2.3.merged.config.default.tsv.gz
    │       │       ├── HEPG2.3.merged.config.outlierZscore.tsv.gz
    │       │       ├── default
    │       │       │   ├── HEPG2.1.barcode_assigned_counts.tsv.gz
    │       │       │   ├── HEPG2.1.barcodesRemoved_assigned_counts.tsv.gz
    │       │       │   ├── HEPG2.1.merged_assigned_counts.tsv.gz
    │       │       │   ├── HEPG2.2.barcode_assigned_counts.tsv.gz
    │       │       │   ├── HEPG2.2.barcodesRemoved_assigned_counts.tsv.gz
    │       │       │   ├── HEPG2.2.merged_assigned_counts.tsv.gz
    │       │       │   ├── HEPG2.3.barcode_assigned_counts.tsv.gz
    │       │       │   ├── HEPG2.3.barcodesRemoved_assigned_counts.tsv.gz
    │       │       │   ├── HEPG2.3.merged_assigned_counts.tsv.gz
    │       │       │   ├── HEPG2.allreps.merged.combined.tsv.gz
    │       │       │   ├── HEPG2.allreps.merged.tsv.gz
    │       │       │   ├── HEPG2.allreps.merged_barcode_assigned_counts.tsv.gz
    │       │       │   ├── HEPG2.allreps_minThreshold.merged.combined.tsv.gz
    │       │       │   ├── HEPG2.allreps_minThreshold.merged.tsv.gz
    │       │       │   └── HEPG2.allreps_minThreshold.merged_barcode_assigned_counts.tsv.gz
    │       │       └── outlierZscore
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
    │       │   └── fromFile.tsv.gz
    │       ├── counts
    │       │   ├── HEPG2.1.DNA.filtered_counts.tsv.gz
    │       │   ├── HEPG2.1.DNA.final_counts.tsv.gz
    │       │   ├── HEPG2.1.RNA.filtered_counts.tsv.gz
    │       │   ├── HEPG2.1.RNA.final_counts.tsv.gz
    │       │   ├── HEPG2.1.merged.config.default.tsv.gz
    │       │   ├── HEPG2.1.merged.config.outlierZscore.tsv.gz
    │       │   ├── HEPG2.2.DNA.filtered_counts.tsv.gz
    │       │   ├── HEPG2.2.DNA.final_counts.tsv.gz
    │       │   ├── HEPG2.2.RNA.filtered_counts.tsv.gz
    │       │   ├── HEPG2.2.RNA.final_counts.tsv.gz
    │       │   ├── HEPG2.2.merged.config.default.tsv.gz
    │       │   ├── HEPG2.2.merged.config.outlierZscore.tsv.gz
    │       │   ├── HEPG2.3.DNA.filtered_counts.tsv.gz
    │       │   ├── HEPG2.3.DNA.final_counts.tsv.gz
    │       │   ├── HEPG2.3.RNA.filtered_counts.tsv.gz
    │       │   ├── HEPG2.3.RNA.final_counts.tsv.gz
    │       │   ├── HEPG2.3.merged.config.default.tsv.gz
    │       │   ├── HEPG2.3.merged.config.outlierZscore.tsv.gz
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
    │       ├── qc_metrics.HEPG2.fromFile.default.json
    │       ├── qc_metrics.HEPG2.fromFile.outlierZscore.json
    │       ├── qc_report.HEPG2.fromFile.default.html
    │       ├── qc_report.HEPG2.fromFile.outlierZscore.html
    │       ├── reporter_experiment.barcode.HEPG2.fromFile.default.all.tsv.gz
    │       ├── reporter_experiment.barcode.HEPG2.fromFile.default.min_oligo_threshold_10.tsv.gz
    │       ├── reporter_experiment.barcode.HEPG2.fromFile.outlierZscore.all.tsv.gz
    │       ├── reporter_experiment.barcode.HEPG2.fromFile.outlierZscore.min_oligo_threshold_10.tsv.gz
    │       ├── reporter_experiment.oligo.HEPG2.fromFile.default.all.tsv.gz
    │       ├── reporter_experiment.oligo.HEPG2.fromFile.default.min_oligo_threshold_10.tsv.gz
    │       ├── reporter_experiment.oligo.HEPG2.fromFile.outlierZscore.all.tsv.gz
    │       ├── reporter_experiment.oligo.HEPG2.fromFile.outlierZscore.min_oligo_threshold_10.tsv.gz
    │       └── statistic
    │           ├── assigned_counts
    │           │   └── fromFile
    │           │       ├── default
    │           │       │   ├── HEPG2.DNA.pairwise.minThreshold.png
    │           │       │   ├── HEPG2.DNA.pairwise.png
    │           │       │   ├── HEPG2.RNA.pairwise.minThreshold.png
    │           │       │   ├── HEPG2.RNA.pairwise.png
    │           │       │   ├── HEPG2.Ratio.pairwise.minThreshold.png
    │           │       │   ├── HEPG2.Ratio.pairwise.png
    │           │       │   ├── HEPG2.average_allreps.merged.tsv.gz
    │           │       │   ├── HEPG2.barcodesPerInsert.png
    │           │       │   ├── HEPG2.dna_vs_rna.png
    │           │       │   ├── HEPG2.dna_vs_rna_minThreshold.png
    │           │       │   ├── HEPG2.group_barcodesPerInsert_box.png
    │           │       │   └── HEPG2.group_barcodesPerInsert_box_minThreshold.png
    │           │       └── outlierZscore
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
    │           │   │   └── fromFile
    │           │   │       ├── HEPG2.default.DNA.perBarcode.png
    │           │   │       ├── HEPG2.default.RNA.perBarcode.png
    │           │   │       ├── HEPG2.default.barcode.DNA.pairwise.png
    │           │   │       ├── HEPG2.default.barcode.RNA.pairwise.png
    │           │   │       ├── HEPG2.default.barcode.Ratio.pairwise.png
    │           │   │       ├── HEPG2.outlierZscore.DNA.perBarcode.png
    │           │   │       ├── HEPG2.outlierZscore.RNA.perBarcode.png
    │           │   │       ├── HEPG2.outlierZscore.barcode.DNA.pairwise.png
    │           │   │       ├── HEPG2.outlierZscore.barcode.RNA.pairwise.png
    │           │   │       └── HEPG2.outlierZscore.barcode.Ratio.pairwise.png
    │           │   └── counts
    │           │       ├── HEPG2.default.DNA.perBarcode.png
    │           │       ├── HEPG2.default.RNA.perBarcode.png
    │           │       ├── HEPG2.default.barcode.DNA.pairwise.png
    │           │       ├── HEPG2.default.barcode.RNA.pairwise.png
    │           │       ├── HEPG2.default.barcode.Ratio.pairwise.png
    │           │       ├── HEPG2.outlierZscore.DNA.perBarcode.png
    │           │       ├── HEPG2.outlierZscore.RNA.perBarcode.png
    │           │       ├── HEPG2.outlierZscore.barcode.DNA.pairwise.png
    │           │       ├── HEPG2.outlierZscore.barcode.RNA.pairwise.png
    │           │       └── HEPG2.outlierZscore.barcode.Ratio.pairwise.png
    │           ├── bc_overlap.assigned_counts.default.fromFile.tsv
    │           ├── bc_overlap.assigned_counts.outlierZscore.fromFile.tsv
    │           ├── bc_overlap.counts.default.tsv
    │           ├── bc_overlap.counts.outlierZscore.tsv
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
    │           ├── statistic_assigned_bc_correlation_merged.fromFile.default.tsv
    │           ├── statistic_assigned_bc_correlation_merged.fromFile.outlierZscore.tsv
    │           ├── statistic_assigned_counts_merged.fromFile.default.tsv
    │           ├── statistic_assigned_counts_merged.fromFile.outlierZscore.tsv
    │           ├── statistic_assigned_counts_single.fromFile.default.tsv
    │           ├── statistic_assigned_counts_single.fromFile.outlierZscore.tsv
    │           ├── statistic_bc_correlation_merged.default.tsv
    │           ├── statistic_bc_correlation_merged.outlierZscore.tsv
    │           ├── statistic_oligo_correlation_merged.fromFile.default.tsv
    │           └── statistic_oligo_correlation_merged.fromFile.outlierZscore.tsv
    └── logs
