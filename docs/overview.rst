.. _Overview:

=====================
Overview
=====================

This pipeline processes sequencing data from Massively Parallel Reporter Assays (MPRAs) to create count tables for candidate sequences tested in the experiment.

This package contains three utilities:

ASSOCIATION
-----------
This utility takes in library association sequencing data (FASTQ) and a design file (FASTA) to assign barcodes to the corresponding elements tested. Functionality includes filtering for quality and coverage of barcodes. This utility must be run before the COUNT utility. See :ref:`Association` for more details.

COUNT 
-----
This utility processes sequence data (FASTQ) of barcodes from the DNA and RNA fractions of the MPRA experiment and outputs count tables labeled with the element tested and a label provided in the design file. This utility can process multiple replicates and conditions in a parallelized manner. Based on a user specified flag, the pipeline will either output normalized activity for each tested sequence, or will combine the results into a single count matrix compatible with MPRAnalyze. See :ref:`Count` for more details.

ASSIGNMENT 
-----
Not documented yet.
