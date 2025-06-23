.. _Overview:

=====================
Overview
=====================

MPRAsnakeflow is a pipeline designed to process sequencing data from Massively Parallel Reporter Assays (MPRAs) to generate count tables for candidate sequences tested in the experiment.

This package contains two main utilities:

ASSIGNMENT (Library Association)
---------------------------------
This utility processes library association sequencing data (FASTQ) and a design file (FASTA) to assign barcodes to the corresponding elements tested. Key features include:
- Assigning barcodes to candidate sequences.
- Filtering for quality and coverage of barcodes.

This utility must be run before the :ref:`Experiment` utility. For more details, see :ref:`Assignment`.

EXPERIMENT (Count Processing)
-----------------------------
This utility processes sequence data (FASTQ) of barcodes from the DNA and RNA fractions of the MPRA experiment and outputs count tables labeled with the tested elements and their corresponding labels from the design file. Key features include:
- Processing multiple replicates and conditions in a parallelized manner.
- Generating normalized activity for each tested sequence or combining results into a single count matrix compatible with MPRAnalyze.

For more details, see :ref:`Experiment`.
