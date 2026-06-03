.. _GSE316891_example:

.. role:: bash(code)
    :language: bash

================================
GSE316891 (Yan et al. L1a1 MPRA)
================================

This example demonstrates running the assignment and experiment workflows on the L1a1 MPRA dataset from Yan et al. Data is available on `GEO:GSE316891 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316891>`_. No preprint or publication is available at the moment, but the data is already public and can be used for testing MPRAsnakeflow.

Prerequisites
=============

This example depends on the following data and software:

Installation of MPRAsnakeflow
-----------------------------
Please install conda, set up the MPRAsnakeflow environment, and clone the current MPRAsnakeflow master branch. See :ref:`Installation` for details.

Assignment Workflow
===================

Download and Prepare Data
-------------------------

.. code-block:: bash

	 mkdir -p data/assignment
	 cd data/assignment
	 curl https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9462nnn/GSM9462019/suppl/GSM9462019%5FL1a1%2Efasta%2Egz | zcat > GSM9462019_L1a1.fa
	 cat GSM9462019_L1a1.fa | sed 's/\[/|/g' | sed 's/\]/|/g' > GSM9462019_L1a1.fix.fa
	 prefetch --max-size 20GB SRR36893758
	 fastq-dump --gzip --split-files SRR36893758

The folder should look like this:

.. code-block:: text

	 data/assignment/
	 ├── GSM9462019_L1a1.fa
	 ├── GSM9462019_L1a1.fix.fa
	 ├── SRR36893758_1.fastq.gz
	 ├── SRR36893758_2.fastq.gz
	 └── SRR36893758_3.fastq.gz

Configure Assignment Workflow
-----------------------------

Create a config file (e.g. ``config_assignment.yaml``) with the following content:

.. code-block:: yaml

	version: "0.7.0"
	assignments:
	  GSE316891L1a1:
	    bc_length: 15
        alignment_tool:
          split_number: 30
          tool: bbmap
        FWD:
          - data/assignment/SRR36893758_1.fastq.gz
        REV:
          - data/assignment/SRR36893758_2.fastq.gz
        BC:
          - data/assignment/SRR36893758_3.fastq.gz
        design_check:
          sequence_start: 1
          sequence_length: 270
        design_file: data/assignment/GSM9462019_L1a1.fix.fa
        configs:
          default: {}

Run Assignment Workflow
-----------------------

.. code-block:: bash

    snakemake -c all --sdm apptainer conda \
    --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile \
    --configfile config_assignment.yaml

Experiment Workflow
===================

Download and Prepare Data
-------------------------

.. code-block:: bash

    mkdir -p data/experiment
    cd data/

    wget -O GSM9462019_L1a1_oligos_to_barcodes.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9462nnn/GSM9462019/suppl/GSM9462019%5FL1a1%5Foligos%5Fto%5Fbarcodes%2Etxt%2Egz
    zcat GSM9462019_L1a1_oligos_to_barcodes.txt.gz | tail -n +2 | awk -F',' -v "OFS=\t" '{print $2,$1}' | sed 's/\[/|/g' | sed 's/\]/|/g' | sort | gzip -c > GSM9462019_L1a1_oligos_to_barcodes.tsv.gz
    cd experiment
    prefetch --max-size 30GB SRR36893764 SRR36893763 SRR36893762 SRR36893761 SRR36893760 SRR36893759
    fastq-dump --gzip --split-files SRR36893764
    fastq-dump --gzip --split-files SRR36893763
    fastq-dump --gzip --split-files SRR36893762
    fastq-dump --gzip --split-files SRR36893761
    fastq-dump --gzip --split-files SRR36893760
    fastq-dump --gzip --split-files SRR36893759

The folder should look like this:

.. code-block:: text

    data/experiment/
    ├── SRR36893759_1.fastq.gz
    ├── SRR36893759_2.fastq.gz
    ├── SRR36893759_3.fastq.gz
    ├── SRR36893760_1.fastq.gz
    ├── SRR36893760_2.fastq.gz
    ├── SRR36893760_3.fastq.gz
    ├── SRR36893761_1.fastq.gz
    ├── SRR36893761_2.fastq.gz
    ├── SRR36893761_3.fastq.gz
    ├── SRR36893762_1.fastq.gz
    ├── SRR36893762_2.fastq.gz
    ├── SRR36893762_3.fastq.gz
    ├── SRR36893763_1.fastq.gz
    ├── SRR36893763_2.fastq.gz
    ├── SRR36893763_3.fastq.gz
    ├── SRR36893764_1.fastq.gz
    ├── SRR36893764_2.fastq.gz
    └── SRR36893764_3.fastq.gz

Create experiment file (``experiments.csv``):

.. code-block:: text

	 Condition,Replicate,DNA_BC_F,DNA_BC_R,DNA_UMI,RNA_BC_F,RNA_BC_R,RNA_UMI
	 L1a1,1,SRR36893764_1.fastq.gz,SRR36893764_2.fastq.gz,SRR36893764_3.fastq.gz,SRR36893763_1.fastq.gz,SRR36893763_2.fastq.gz,SRR36893763_3.fastq.gz
	 L1a1,2,SRR36893762_1.fastq.gz,SRR36893762_2.fastq.gz,SRR36893762_3.fastq.gz,SRR36893761_1.fastq.gz,SRR36893761_2.fastq.gz,SRR36893761_3.fastq.gz
	 L1a1,3,SRR36893760_1.fastq.gz,SRR36893760_2.fastq.gz,SRR36893760_3.fastq.gz,SRR36893759_1.fastq.gz,SRR36893759_2.fastq.gz,SRR36893759_3.fastq.gz

Configure Experiment Workflow
-----------------------------

Create a config file (e.g. ``config_experiment.yaml``):

.. code-block:: yaml

    version: "0.7"
    experiments:
      GSE316891L1a1:
        bc_length: 15
        umi_length: 16
        data_folder: data/experiment
        experiment_file: experiments.csv
        assignments:
          L1a1Assignment:
            type: file
            assignment_file: results/assignment/GSE316891L1a1/assignment_barcodes.default.tsv.gz
          GEOAssignment:
            type: file
            assignment_file: data/GSM9462019_L1a1_oligos_to_barcodes.tsv.gz
        configs:
          default: {}

Run Experiment Workflow
-----------------------

.. code-block:: bash

    snakemake -c all --sdm apptainer conda \
    --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile \
    --configfile config_experiment.yaml

Results
=======

Assignment output files will be in :code:`results/assignment/GSE316891L1a1`. The final assignment is :code:`results/assignment/GSE316891L1a1/assignment_barcodes.default.tsv.gz` and the QC report is :code:`results/assignment/GSE316891L1a1/qc_report.default.html`. We also provide the `Assignment QC report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE316891L1a1.assignment.qc_report.default.html>`_.

Inspecting their assignment file, they retrieve 77795 out of 79812 oligos with 113 mean and 93 median barcodes per oligo. Our default assignment retrieves 77208 oligos with 87 mean and 71 median barcodes per oligo. We are very strict with the majority vote (75% of the reads need to support the same oligo), which leads to a lower number of assigned oligos, but also a higher confidence in the assignment. You can adjust the config to be less strict and retrieve more oligos, but we recommend inspecting the QC report and assignment statistics before doing so.


Experiment output files will be in :code:`results/experiments/GSE316891L1a1`. The main count files are :code:`results/experiments/GSE316891L1a1/reporter_experiment.barcode.L1a1.GEOAssignment.default.all.tsv.gz` and :code:`results/experiments/GSE316891L1a1/reporter_experiment.barcode.L1a1.L1a1Assignment.default.all.tsv.gz`. QC reports are in the same folder and can be inspected for our `L1a1Assignment here <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE316891.experiment.qc_report.L1a1Assignment.default.html>`_ and for their `GEOAssignment here <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE316891.experiment.qc_report.GEOAssignment.default.html>`_. QC reports are very similar for both assignments with slightly more barcodes assigned in the GEOAssignment due to the richer assignment file. The final count files are very similar between the two assignments. Here we report the correlations of the normalized DNA, normalized RNA, and the log2 fold change between RNA and DNA of the two assignments, separated by replicate:


.. list-table:: Correlations Between L1a1Assignment and GEOAssignment
   :header-rows: 1

   * - replicate
     - number of shared oligos
     - pearson log2FC
     - spearman log2FC
     - pearson RNA
     - spearman RNA
     - pearson DNA
     - spearman DNA
   * - 1
     - 75388
     - 0.88
     - 0.89
     - 0.89
     - 0.92
     - 0.82
     - 0.83
   * - 2
     - 75462
     - 0.88
     - 0.9
     - 0.91
     - 0.93
     - 0.81
     - 0.83
   * - 3
     - 75151
     - 0.87
     - 0.88
     - 0.9
     - 0.91
     - 0.82
     - 0.84
