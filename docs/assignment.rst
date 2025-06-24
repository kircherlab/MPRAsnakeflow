.. _Assignment:

=====================
Assignment
=====================

.. image:: MPRAsnakeflow_assignment_overview.png

Input Files
===============

Fastq Files
-----------
- 2-3 Fastq files from library association sequencing
- Candidate regulatory sequence (CRS) sequencing, forward and reverse read (paired-end)
- (Optional) Index read with barcode (see read structure figure, Panel A). The barcode (BC) can also be present at the beginning of the forward read followed by a linker (Panel B).
- (Optional) Adapter trimming from the single end read that should contains the sequenced oligo (provided by a single forward read or merged between forward and reverse reads). It is highly recommended that final processed read has the same length as the designed oligos and contain only oligo information.

Design File
-----------
Multi-FASTA file of CRS sequences with unique headers describing each tested sequence.

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

.. note:: Headers of the design file must be unique (before any space), cannot contain :code:`[` or :code:`]`, and should not contain duplicated sequences (forward and antisense).

Config File
-----------
Multiple mapping strategies are implemented to find the corresponding CRS sequence for each read. The mapping strategy can be chosen in the config file (bbmap, bwa mem, or exact matches). The config file also defines the filtering of the mapping results and is a YAML file.

Example of an assignment file using bbmap and the standard filtering (we recommend using bbmap as the default):

.. literalinclude:: ../config/example_assignment_bbmap.yaml
   :language: yaml

Example of an assignment file using bwa and the standard filtering:

.. literalinclude:: ../config/example_assignment_bwa.yaml
   :language: yaml

Example of an assignment file using exact matches with non-default filtering of barcodes:

.. literalinclude:: ../config/example_assignment_exact_lazy.yaml
   :language: yaml

Example of an assignment file using exact matches and read 1 with BC, linker, and oligo (no separate BC index read):

.. literalinclude:: ../config/example_assignment_exact_linker.yaml
   :language: yaml

If you want to use the strand sensitivity option (e.g., testing enhancers in both directions), you can add the following to the config file: :code:`strand_sensitive: {enable: true}`. Otherwise, MPRAsnakeflow will give you an error because it cannot handle the same sequences in both sense and antisense directions. This is an issue with the mappers because they do not consider the strand and will always call your read ambiguous due to multiple matches.

Snakemake
============================

Options
---------------

With :code:`--help` or :code:`-h`, you can see the help message.

Mandatory arguments:
  :\-\-cores:                 
    Use at most N CPU cores/jobs in parallel. If N is omitted or 'all', the limit is set to the number of available CPU cores. In the case of cluster/cloud execution, this argument sets the number of total cores used over all jobs (made available to rules via workflow.cores). (default: None)
  :\-\-configfile:
    Specify or overwrite the config file of the workflow (see the docs). Values specified in JSON or YAML format are available in the global config dictionary inside the workflow. Multiple files overwrite each other in the given order. Missing keys in previous config files are extended by following config files. (default: None)
  :\-\-sdm:             
    **Required to run MPRAsnakeflow.** : :code:`--sdm conda` or :code:`--sdm apptainer conda`. Uses the defined conda environment per rule. We highly recommend using apptainer, where we build a predefined Docker container with all software installed within it. :code:`--sdm conda` installs the conda environments during the first execution of the workflow. If this flag is not set, the conda/apptainer directive is ignored. (default: False)

Recommended arguments:
  :\-\-snakefile:             
    You should not need to specify this. By default, Snakemake will search for 'Snakefile', 'snakefile', 'workflow/Snakefile', or 'workflow/snakefile' beneath the current working directory, in this order. Only if you definitely want a different layout, you need to use this parameter. This is very useful when you want to have the results in a different folder than MPRAsnakeflow is in. (default: None)

Useful arguments:
  :-n:                      
    Do not execute anything, and display what would be done. If you have a very large workflow, use --dry-run --quiet to just print a summary of the DAG of jobs. (default: False)
  :\-\-touch, -t:             
    Touch output files (mark them up to date without really changing them) instead of running their commands. This is used to pretend that the rules were executed, in order to fool future invocations of Snakemake. Fails if a file does not yet exist. Note that this will only touch files that would otherwise be recreated by Snakemake (e.g., because their input files are newer). For enforcing a touch, combine this with --force, --forceall, or --forcerun. Note, however, that you lose the provenance information when the files have been created in reality. Hence, this should be used only as a last resort. (default: False)

Rules
-------------

Rules run by Snakemake in the assignment utility:

- **all**: The overall rule. Defines what final output files are expected.
- **all_assignments**: Run all steps of the assignment workflow.
- **assignment_attach_idx**: Extract the index sequence and add it to the header.
- **assignment_collect**: Collect mapped reads into one BAM.
- **assignment_collectBCs**: Get the barcodes.
- **assignment_fastq_split**: Split the fastq files into N files for parallelization. N is given by `split_read` in the configuration file.
- **assignment_filter**: Filter the barcodes file based on the config given in the config file. Results are here: :code:`results/assignment/<assignment_name>/assignment_barcodes.<config_name>.tsv.gz`.
- **assignment_flagstat**: Run samtools flagstat. Results are in :code:`results/assignment/<assignment_name>/statistic/assignment/bam_stats.txt`.
- **assignment_mapping_bwa**: Map the reads to the reference using BWA.
- **assignment_merge**: Merge the forward, reverse, and barcode fastq files into one.

Output
==========

The output can be found in the folder defined by the option :code:`results/assignment/`. It is structured in folders of the condition as follows:

Files
-------------
Once the pipeline is finished running, all the output files can be seen in the results folder. This pipeline also generates a QC report. 

For more details, refer to the `HTML QC report <https://kircherlab.github.io/mprasnakeflow/assignment.html>`_.

File tree of the result folder (names in :code:`< >` can be specified in the config file):

.. code-block:: text

    ├── assignment
    └── <assignment_name>
        ├── BCs
        ├── aligned_merged_reads.bam
        ├── aligned_merged_reads.bam.bai
        ├── assignment_barcodes.default.tsv.gz
        ├── assignment_barcodes_with_ambiguous.default.tsv.gz
        ├── barcodes_incl_other.tsv.gz
        ├── bbmap
        ├── design_check.done
        ├── design_check.err
        ├── fastq
        │   └── splits
        ├── qc_report.default.html
        ├── reference
        │   └── reference.fa
        └── statistic
            ├── assigned_counts.default.tsv
            ├── assignment
            │   └── bam_stats.txt
            ├── assignment.default.png
            ├── assignment.default.tsv.gz
            └── total_counts.tsv

Key output files:
- **qc_report.default.html**: QC report of the assignment.
- **total_counts.tsv**: Raw statistics of barcodes mapped to oligos.
- **assigned_counts.<config_name>.tsv**: Statistics of barcodes mapped to oligos after filtering.
- **assignment.<config_name>.tsv.gz**: Average/median support of barcodes per oligo.
- **reference.fa**: Design file.
- **aligned_merged_reads.bam**: Sorted BAM file for oligo alignment.
- **assignment_barcodes.<config_name>.tsv.gz**: Mapping file of barcodes to sequences.
