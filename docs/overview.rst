.. _Overview:

=====================
Overview
=====================

This pipeline processes sequencing data from Massively Parallel Reporter Assays (MPRAs) to create count tables for candidate sequences tested in the experiment.

This package contains three utilities:

ASSIGNMENT (previous associations)
----------------------------------
This utility takes in library association sequencing data (FASTQ) and a design file (FASTA) to assign barcodes to the corresponding elements tested. Functionality includes filtering for quality and coverage of barcodes. This utility must be run before the :ref:`Experiment` utility. See :ref:`Assignment` for more details.

EXPERIMENT (previous count)
---------------------------
This utility processes sequence data (FASTQ) of barcodes from the DNA and RNA fractions of the MPRA experiment and outputs count tables labeled with the element tested and a label provided in the design file. This utility can process multiple replicates and conditions in a parallelized manner. Based on a user specified flag, the pipeline will either output normalized activity for each tested sequence, or will combine the results into a single count matrix compatible with MPRAnalyze. See :ref:`Experiment` for more details.
