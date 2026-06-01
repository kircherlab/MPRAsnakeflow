.. _Examples overview:

=================
Examples Overview
=================

This page summarizes all example datasets in MPRAsnakeflow and links to the full step-by-step instructions for each workflow.

Core Workflow Examples
======================

* :doc:`assignment_example1`

  Basic assignment-only workflow on 5'/5' WT MPRA data in HepG2 from `Klein et al. (2019) <http://dx.doi.org/10.1038/s41592-020-0965-y>`_. Use this example to learn barcode-to-oligo assignment from raw reads.

* :doc:`count_example1`

  Basic experiment/count workflow on 5'/5' WT MPRA data in HepG2 from `Klein et al. (2019) <http://dx.doi.org/10.1038/s41592-020-0965-y>`_, using a precomputed assignment file.

* :doc:`combined_example1`

  Combined assignment + experiment workflow on the same HepG2 WT MPRA dataset from `Klein et al. (2019) <http://dx.doi.org/10.1038/s41592-020-0965-y>`_, useful as an end-to-end minimal example.

Published Dataset Examples
==========================

* :doc:`plasmid_example`

  ENCODE plasmid-based MPRA in A549 from Tewhey lab, published in `Gosai et al. (Nature 2024) <https://doi.org/10.1038/s41586-024-08070-z>`_. Demonstrates assignment preprocessing when barcodes are attached to forward reads and experiment setup with shared DNA input across replicates.

* :doc:`complex_readstructure_example`

  Complex-read-structure MPRA from `Abell et al. (Science 2022) <https://doi.org/10.1126/science.abj5117>`_, GEO `GSE174534 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174534>`_. Focuses on handling non-trivial read layouts, trimming/adapters, strand-sensitive assignment, and reverse-complement design handling.

* :doc:`GSE306816_example`

  Deep perturbation STR MPRA from `Zhang et al. (bioRxiv 2025) <https://doi.org/10.1101/2025.09.14.676153>`_, GEO `GSE306816 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306816>`_. Shows assignment and counting for repeat-rich constructs and comparison against published assignment files.

* :doc:`GSE293036_example`

  Multiple sclerosis variant MPRA from `Granitto et al. (G3 2025) <https://doi.org/10.1093/g3journal/jkaf192>`_, GEO `GSE293036 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293036>`_. Demonstrates conversion from supplementary design tables, assignment generation, and condition-specific experiment counting.

* :doc:`GSE316891_example`

  L1a1 MPRA dataset from Yan et al., GEO `GSE316891 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316891>`_. A peer-reviewed paper or preprint link is not currently available in the example source data. Includes assignment generation, experiment counting, and direct comparison between workflow-generated and GEO-provided assignment files.

* :doc:`GSE284330_example`

  Processed HepG2 sub1 MPRA data from `Zaratiana et al. <https://www.cell.com/molecular-cell/fulltext/S1097-2765(26)00232-7>`_, GEO `GSE284330 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284330>`_. STARR-seq like assay, therefore not optimal for the workflow. But it uses designed oligonucleotides we demonstrate experiment-only processing using oligonucleotides as barcodes and shared DNA input across RNA replicates.

* :doc:`GSE325670_example`

  Variation/saturation mutagenesis MPRA from `Hauser et al. (2026) <https://doi.org/10.64898/2026.03.06.710116>`_, GEO umbrella `GSE325670 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325670>`_ (example uses GSE325256). Highlights challenging library assignment settings, including bbmap/bwa-additional-filtering, strand-sensitive assignment, and downstream comparison in experiment counting.
