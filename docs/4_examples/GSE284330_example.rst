.. _GSE284330 example:

.. role:: bash(code)
   :language: bash

=================================================
GSE284330 (Zaratiana et al. STARR-seq like assay)
=================================================

This example runs the experiment/count workflow on processed data from `Zaratiana et al. <https://www.cell.com/molecular-cell/fulltext/S1097-2765(26)00232-7>`_. The data were published in `GEO:GSE284330 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284330>`_.

The GEO series describes a genome-wide MPRA in mouse liver and human HepG2 cells to identify functional cis-regulatory elements and their response to gut microbiota signals. Here we use the processed HepG2 sub1 example.

.. note:: STARR-seq like assay, therefore not optimal for the workflow. But it uses designed oligonucleotides we demonstrate experiment-only processing using oligonucleotides as barcodes (length of 200).

Prerequisites
=============

This example depends on the following data and software:

Installation of MPRAsnakeflow
-----------------------------
Please install conda, set up the MPRAsnakeflow environment, and clone the current MPRAsnakeflow master branch. See :ref:`Installation` for details.

Processed Assignment and Count Data
===================================

This is an experiment-only example. We use a published barcode mapping together with raw count FASTQ files from GEO.

Download and Prepare Data
-------------------------

The example uses HepG2 sub1 DNA and RNA data. The barcode mapping is created from Supplementary Table 1 of the publication.

.. code-block:: bash

    mkdir -p data/experiment
    cd data/experiment

    # Save Supplementary Table 1 as design.tsv, then create the barcode mapping.
    awk -F'\t' -v "OFS=_" '{ $6=toupper($6); print $6"\t"$1,$2,$3,$4,$5}' design.tsv | tail -n +2 | gzip -c > bc_mapping.tsv.gz

    prefetch --max-size 30GB SRR31723721 SRR31723722 SRR31723723 SRR31723720 SRR31723719
    fastq-dump --gzip --split-files SRR31723721 && \
    fastq-dump --gzip --split-files SRR31723722 && \
    fastq-dump --gzip --split-files SRR31723723 && \
    fastq-dump --gzip --split-files SRR31723720 && \
    fastq-dump --gzip --split-files SRR31723719

The folder should look like this:

.. code-block:: text

    data/
    └── experiment
        ├── SRR31723719_1.fastq.gz
        ├── SRR31723719_2.fastq.gz
        ├── SRR31723720_1.fastq.gz
        ├── SRR31723720_2.fastq.gz
        ├── SRR31723721_1.fastq.gz
        ├── SRR31723721_2.fastq.gz
        ├── SRR31723722_1.fastq.gz
        ├── SRR31723722_2.fastq.gz
        ├── SRR31723723_1.fastq.gz
        └── SRR31723723_2.fastq.gz

Reads count data
----------------

.. list-table:: HepG2 sub1 data
   :header-rows: 1

   * - Sample type
     - Replicate
     - Runs
   * - DNA
     - shared input
     - SRR31723721, SRR31723722, SRR31723723
   * - RNA
     - 1
     - SRR31723720
   * - RNA
     - 2
     - SRR31723719

MPRAsnakeflow
=============

We run the experiment/count workflow only, using the processed barcode mapping derived from the publication supplementary table.

Create config files
-------------------

Create the :code:`experiments.csv` file to map DNA/RNA counts to the correct replicates:

.. code-block:: bash

    cat << 'EOF' > experiments.csv
    Condition,Replicate,DNA_BC_F,DNA_BC_R,RNA_BC_F,RNA_BC_R
    HepG2Sub1,1,SRR31723721_1.fastq.gz;SRR31723722_1.fastq.gz;SRR31723723_1.fastq.gz,SRR31723721_2.fastq.gz;SRR31723722_2.fastq.gz;SRR31723723_2.fastq.gz,SRR31723720_1.fastq.gz,SRR31723720_2.fastq.gz
    HepG2Sub1,2,SRR31723721_1.fastq.gz;SRR31723722_1.fastq.gz;SRR31723723_1.fastq.gz,SRR31723721_2.fastq.gz;SRR31723722_2.fastq.gz;SRR31723723_2.fastq.gz,SRR31723719_1.fastq.gz,SRR31723719_2.fastq.gz
    EOF


Create a config file (e.g. ``config_experiment.yaml``). Please note that we have to set the barcode length to 200, which is the length of the designed oligonucleotides used as barcodes in this assay. Because the implemented paird-end BC merging is not optimized for such long reads we use NGmerge ``merge_tool: NGmerge`` and extend the maximum overlap. Also note that we have to use a barcode threshold of 1 because there are not really barcodes. We do a one to one mapping of the oligos:

.. code-block:: yaml

    ---
    version: "0.7"
    experiments:
      GSE284330:
        bc_length: 200
        merge_tool: NGmerge
        NGmerge:
          min_overlap: 80
        adapters:
          FWD:
            five_prime:
              - ATGAG
          REV:
            five_prime:
              - TCGTG
        data_folder: data/experiment
        experiment_file: experiments.csv
        assignments:
          GSE284330:
            type: file
            assignment_file: data/experiment/bc_mapping.tsv.gz
        configs:
          BCthreshold1:
            filter:
              bc_threshold: 1

Run snakemake
-------------

Now we are ready to run MPRAsnakeflow. We run this example on one node with 30 cores.

First, do a dry run with snakemake using :code:`-n`:

.. code-block:: bash

    conda activate mprasnakeflow
    snakemake -c 1 --sdm apptainer conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile config_experiment.yaml -n --quiet rules

If the dry run finishes without errors, run the workflow:

.. code-block:: bash

    snakemake -c 30 --sdm apptainer conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile config_experiment.yaml

.. note:: Please adapt your setup when running in a cluster environment. An example SLURM profile is available in MPRAsnakeflow under :code:`profiles/default/config.yaml`. You can use it with snakemake via :code:`--workflow-profile $PIPELINE/profiles/default`, but adjust it first, especially the :code:`slurm_partition` setting.

Results
=======

All output files are written to :code:`results/experiments/GSE284330`.

The final count file is:

* :code:`results/experiments/GSE284330/reporter_experiment.barcode.HepG2Sub1.GSE284330.BCthreshold1.all.tsv.gz`

You should also inspect the QC report in the same directory, or have a look `here <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE284330.experiment.qc_report.HepG2Sub1.BCthreshold1.html>`_. As expected we have a meadian BC of 1. We also see that we only get a library coverage of less then 20%. This is expected because this assay is not really designed for barcode-based counting. We are doing a exact mapping of 200bp sequenced oligos. The assay is more similar to a STARR-seq like assay where the oligos themselves are sequenced and counted. Therefore, we have a lot of oligos that are not covered by any read and many reads that do not match any oligo. This example is therefore not ideal for MPRAsnakeflow, but it demonstrates how to use the workflow can be used in that satting and highlits potential optimizations for future versions of the workflow, like implementing a mapping strategy for the experiment sub-workflow.
