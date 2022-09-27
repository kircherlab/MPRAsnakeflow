.. _Assignment:

=====================
Assignment
=====================

.. image:: Association_util.png

Input files
===============

Fastq Files
-----------
- 2-3 Fastq files from library association sequencing
- Candidate regulatory sequence (CRS) sequencing, 1 forawrd read and an optional reverse read if paired end sequencing was used
- Barcode sequence, 1 read covering the barcode

Design File
-----------
Fasta file of of CRS sequences with unique headers describing each tested sequence

Example file:

.. code-block:: text

    >CRS1
    GACGGGAACGTTTGAGCGAGATCGAGGATAGGAGGAGCGGA
    >CRS2
    GGGCTCTCTTATATTAAGGGGGTGTGTGAACGCTCGCGATT
    >CRS3
    GGCGCGCTTTTTCGAAGAAACCCGCCGGAGAATATAAGGGA
    >CRS4
    TTAGACCGCCCTTTACCCCGAGAAAACTCAGCTACACACTC


snakemake
============================
 
Options
---------------

With :code:`--help` or :code:`-h` you can see the help message.

Mandatory arguments:
  :\-\-cores:                 
    Use at most N CPU cores/jobs in parallel. If N is omitted or 'all', the limit is set to the number of available CPU cores. In case of cluster/cloud execution, this argument sets the number of total cores used over all jobs (made available to rules via workflow.cores).(default: None)
  :\-\-config, -c:            
    Set or overwrite values in the workflow config object. The workflow config object is accessible as variable config inside the workflow. Default values can be set by providing a JSON file (see Documentation). (default: None)
  :\-\-use-conda:             
    **Required to run MPRAsnakeflow.** If defined in the rule, run job in a conda environment. If this flag is not set, the conda directive is ignored. (default: False)
Recommended arguments:
  :\-\-snakefile:             
    You should not need to specify this. By default, Snakemake will search for 'Snakefile', 'snakefile', 'workflow/Snakefile','workflow/snakefile' beneath the current working directory, in this order. Only if you definitely want a different layout, you need to use this parameter. This is very usefull when you want to have the results in a different folder than MPRAsnakeflow is in. (default: None)
Usefull arguments:
  :-n:                      
    Do not execute anything, and display what would be done. If you have a very large workflow, use --dry-run --quiet to just print a summary of the DAG of jobs. (default: False)
  :\-\-touch, -t:             
    Touch output files (mark them up to date without really changing them) instead of running their commands. This is used to pretend that the rules were executed, in order to fool future invocations of snakemake. Fails if a file does not yet exist. Note that this will only touch files that would otherwise be recreated by Snakemake (e.g. because their input files are newer). For enforcing a touch, combine this with --force, --forceall, or --forcerun. Note however that you loose the provenance information when the files have been created in realitiy. Hence, this should be used only as a last resort. (default: False)

Rules
-------------

Rules run by snakemake in the assignment utility. Some rules will be run only if certain options used and are marked below.

count_bc or count_bc_nolab (if no label file is provided)
  Removes any illegal characters (defined by Piccard) in the label file and design file. Counts the number of reads in the fastq file.

create_BWA_ref
  Creates a BWA reference based on the design file

PE_merge (if paired end fastq files provided)
  Merges the forward and reverse reads covering the CRS using fastq-join

align_BWA_PE or align_BWA_S (if single end mode)
  Uses BWA to align the CRS fastq files to the reference created from the Design File. This will be done for each fastq file chunk based on the split option.

collect_chunks
  merges all bamfiles from each separate alignment

map_element_barcodes
  Assign barcodes to CRS and filters barcodes by user defined parameters for coverage and mapping percentage

filter_barcodes
  Visualize results

.. todo:: These are not the correct files for each cindition in the experiment workflow

Output
==========

The output can be found in the folder defined by the option :code:`results/assignments/`. It is structured in folders of the condition as

Files
-------------

File tree

.. code-block:: text



.. todo:: File tree for the assignment workflow


count_fastq.txt
    number of barcode reads
count_merged.txt
    number of aligned CRS reads
design_rmIllegalChars.fa
    Design file with illegal characters removed
label_rmIllegalChars.txt
    Label file with illegal characters removed
s_merged.bam
    sorted bamfile for CRS alignment
${name}_coords_to_barcodes.pickle
    pickle file containing a python dictionary of CRS/barcode mappings
\*.png
    Visualization of number of barcodes mapping to enhancers

.. todo:: These are not the correct filesin the experiment workflow
