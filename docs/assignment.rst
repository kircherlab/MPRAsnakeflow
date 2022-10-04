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
  :\-\-configfile:
    Specify or overwrite the config file of the workflow (see the docs). Values specified in JSON or YAML format are available in the global config dictionary inside the workflow. Multiple files overwrite each other in the given order. Thereby missing keys in previous config files are extended by following configfiles. Note that this order also includes a config file defined in the workflow definition itself (which will come first). (default: None)
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

Rules run by snakemake in the assignment utility.

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
   Run samtools flagstat. Results are in :code:`results/assignment/<assignment_name>/statistic/assignment/bam_stats.txt`
assignment_getBCs
   Get the barcodes (not filtered). Results are in :code:`results/assignment/<assignment_name>/barcodes_incl_other.sorted.tsv.gz`
assignment_statistic_totalCounts
   Statistic of the total (unfiltered counts). Results are in :code:`results/assignment/<assignment_name>/statistic/total_counts.tsv.gz`
assignment_filter
   Filter the barcodes file based on the config given in the config-file. Results for this run are here :code:`results/assignment/<assignment_name>/assignment_barcodes.<config_name>.sorted.tsv.gz`.
assignment_statistic_assignedCounts
   Statistic of filtered the assigned counts. Result is here :code:`results/assignment/<assignment_name>/statistic/assigned_counts.<config_name>.tsv.gz`.
assignment_statistic_assignment
   Statistic of the filtered assignment.  Result is here :code:`results/assignment/<assignment_name>/statistic/assignment.<config_name>.tsv.gz` and a plot here :code:`results/assignment/<assignment_name>/statistic/assignment.<config_name>.png`.


Output
==========

The output can be found in the folder defined by the option :code:`results/assignment/`. It is structured in folders of the condition as

Files
-------------

File tree of the result folder (names in :code:`< >` can be specified in the config file.)

.. code-block:: text

    ├── assignment
    │   └── <assignment_name>
    │       ├── aligned_merged_reads.bam
    │       ├── aligned_merged_reads.bam.bai
    │       ├── assignment_barcodes.<config_name>.sorted.tsv.gz
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
    │           ├── assigned_counts.<config_name>.tsv.gz
    │           ├── assignment
    │           │   └── bam_stats.txt
    │           ├── assignment.<config_name>.png
    │           ├── assignment.<config_name>.tsv.gz
    │           └── total_counts.tsv.gz




total_counts.tsv.gz
    Statistic of BCs mapped to oligos.
assigned_counts.<config_name>.tsv.gz
    Statistic of BCs mapped to oligos after fitering defined by config.
assignment.<config_name>.tsv.gz
    Average/median support of BC per oligo. Oligos with >= 15 BCs.
reference.fa
    Design file.
aligned_merged_reads.bam
    Sorted bamfile for oligo alignment
barcodes_incl_other.sorted.tsv.gz
    Complete list of all barcodes found in mapping file (ambigous and unambigous) with mappings (if possible)
assignment_barcodes.<config_name>.sorted.tsv.gz
    Mapping file of barcodes to sequence. 
assignment.<config_name>.png
    Visualization of number of barcodes mapping to oligo.
bam_stats.txt
    samtools bamstat output.
