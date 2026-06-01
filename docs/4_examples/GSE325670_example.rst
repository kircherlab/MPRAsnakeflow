.. _GSE325670 example:

.. role:: bash(code)
   :language: bash

=========================================
GSE325670 (Hauser et al. Variant MPRA)
=========================================

This example documents the `GSE325256 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325256>`_ variation and saturation mutagenesis MPRA from `Hauser et al. <https://doi.org/10.64898/2026.03.06.710116>`_, which is one of the three MPRA datasets published under GEO accession `GSE325670 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325670>`_. The other two datasets are `GSE324616 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324616>`_ and `GSE325585 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325585>`_.

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

    # download design file
    curl https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9640nnn/GSM9640449/suppl/GSM9640449%5FVariantMPRAdesignFile%2Efa%2Egz | zcat > GSM9640449_VariantMPRAdesignFile.fa

    # Keep only the first entry for duplicate ID + sequence records in the multi-FASTA.
    awk 'BEGIN { RS=">"; ORS="" }
    NR==1 { next }
    {
        n = split($0, lines, "\n")
        id = lines[1]
        seq = ""
        body = ""
        for (i = 2; i <= n; i++) {
        if (lines[i] == "") continue
        seq = seq lines[i]
        body = (body == "" ? lines[i] : body "\n" lines[i])
        }
        key = id "\t" seq
        if (!(key in seen)) {
        seen[key] = 1
        print ">" id "\n" body "\n"
        }
    }' GSM9640449_VariantMPRAdesignFile.fa > GSM9640449_VariantMPRAdesignFile.dedup.fa

    # Merge identical sequences and trim 15 nt from each end.
    awk 'BEGIN { RS=">"; ORS="" }
    NR==1 { next }
    {
        n = split($0, lines, "\n")
        id = lines[1]
        seq = ""
        for (i = 2; i <= n; i++) {
        if (lines[i] == "") continue
        seq = seq lines[i]
        }
        if (seq == "") next

        if (!(seq in seen)) {
        seen[seq] = 1
        order[++m] = seq
        headers[seq] = id
        } else {
        headers[seq] = headers[seq] "-" id
        }
    }
    END {
        for (i = 1; i <= m; i++) {
        s = order[i]
        if (length(s) <= 30) continue
        s_trim = substr(s, 16, length(s) - 30)
        print ">" headers[s] "\n" s_trim "\n"
        }
    }' GSM9640449_VariantMPRAdesignFile.dedup.fa > GSM9640449_VariantMPRAdesignFile.seq_merged.fa


    # download sequencing data
    wget -O SRR37895574_1.fastq.gz https://sra-pub-src-2.s3.amazonaws.com/SRR37895574/RMHex03_bcassoc_R1_trimmed.fastq.gz.1
    wget -O SRR37895574_2.fastq.gz https://sra-pub-src-2.s3.amazonaws.com/SRR37895574/RMHex03_bcassoc_R2_trimmed.fastq.gz.1
    wget -O SRR37895574_3.fastq.gz https://sra-pub-src-2.s3.amazonaws.com/SRR37895574/RMHex03_bcassoc_R3_trimmed.fastq.gz.1

    cd ../../


The folder should look like this:

.. code-block:: text

    assignment
    ├── GSM9640449_VariantMPRAdesignFile.dedup.fa
    ├── GSM9640449_VariantMPRAdesignFile.fa
    ├── GSM9640449_VariantMPRAdesignFile.seq_merged.fa
    ├── SRR37895574_1.fastq.gz
    ├── SRR37895574_2.fastq.gz
    └── SRR37895574_3.fastq.gz


Configure Assignment Workflow
-----------------------------

Create a config file (e.g. ``config_assignment.yaml``) with the following content:

.. code-block:: yaml

    ---
    version: "0.7.0"
    assignments:
      GSE325256:
        bc_length: 15
        alignment_tool:
          split_number: 30
          tool: bbmap
        FWD:
          - data/assignment/SRR37895574_1.fastq.gz
        BC:
          - data/assignment/SRR37895574_2.fastq.gz
        REV:
          - data/assignment/SRR37895574_3.fastq.gz
        design_file: data/assignment/GSM9640449_VariantMPRAdesignFile.seq_merged.fa
        design_check:
          sequence_collisions: false
        configs:
          default: {}
          lowFraction:
            fraction: 0.51
      GSE325256BWA:
        bc_length: 15
        alignment_tool:
          split_number: 30
          tool: bwa-additional-filtering
          configs:
            min_mapping_quality: 1
            sequence_length:
              min: 200
              max: 225
            alignment_start:
              min: 1
              max: 3
        FWD:
          - data/assignment/SRR37895574_1.fastq.gz
        BC:
          - data/assignment/SRR37895574_2.fastq.gz
        REV:
          - data/assignment/SRR37895574_3.fastq.gz
        design_file: data/assignment/GSM9640449_VariantMPRAdesignFile.seq_merged.fa
        design_check:
          sequence_collisions: false
        configs:
          default: {}
          lowFraction:
            fraction: 0.51
      GSE325256CIGAR:
        bc_length: 15
        alignment_tool:
          split_number: 30
          tool: bbmap
          configs:
            min_mapping_quality: 0
            cigar_filter_regex: 208M|210M|212M|214M|216M|218M|219M
        FWD:
          - data/assignment/SRR37895574_1.fastq.gz
        BC:
          - data/assignment/SRR37895574_2.fastq.gz
        REV:
          - data/assignment/SRR37895574_3.fastq.gz
        design_file: data/assignment/GSM9640449_VariantMPRAdesignFile.seq_merged.fa
        design_check:
          sequence_collisions: false
        configs:
          default: {}
          lowFraction:
            fraction: 0.51
      GSE325256StrandSensitive:
        bc_length: 15
        strand_sensitive:
          enable: true
        alignment_tool:
          split_number: 30
          tool: bbmap
        FWD:
          - data/assignment/SRR37895574_1.fastq.gz
        BC:
          - data/assignment/SRR37895574_2.fastq.gz
        REV:
          - data/assignment/SRR37895574_3.fastq.gz
        design_file: data/assignment/GSM9640449_VariantMPRAdesignFile.seq_merged.fa
        design_check:
          sequence_collisions: false
        configs:
          default: {}
          lowFraction:
            fraction: 0.51

Run Assignment Workflow
-----------------------

.. code-block:: bash

    conda activate mprasnakeflow
    snakemake -c 1 --sdm conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile config_assignment.yaml -n -q --set-threads assignment_mapping_bbmap=10 --resources mem_mb=60000
    snakemake -c 30 --sdm conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile config_assignment.yaml -q --set-threads assignment_mapping_bbmap=10 --resources mem_mb=60000


Experiment Workflow
===================

Download and Prepare Data
-------------------------

.. list-table:: Reads count data
   :header-rows: 1

   * - Replicate
     - DNA SRA ID
     - RNA SRA ID
   * - 1
     - SRR37425536
     - SRR37425544
   * - 2
     - SRR37425537
     - SRR37425545
   * - 3
     - SRR37425531
     - SRR37425546
   * - 4
     - SRR37425532
     - SRR37425533

.. code-block:: bash


    mkdir -p data/experiment
    cd data/experiment
    # download data
    prefetch --max-size 30GB SRR37425536 SRR37425544 SRR37425537 SRR37425545 SRR37425531 SRR37425546 SRR37425532 SRR37425533;
    fastq-dump --gzip --split-files SRR37425536 && \
    fastq-dump --gzip --split-files SRR37425544 && \
    fastq-dump --gzip --split-files SRR37425537 && \
    fastq-dump --gzip --split-files SRR37425545 && \
    fastq-dump --gzip --split-files SRR37425531 && \
    fastq-dump --gzip --split-files SRR37425546 && \
    fastq-dump --gzip --split-files SRR37425532 && \
    fastq-dump --gzip --split-files SRR37425533

    cd ../
    wget -O GSM9640449_VariantMPRAbarcodeDictionary.tsv.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9640nnn/GSM9640449/suppl/GSM9640449%5FVariantMPRAbarcodeDictionary%2Etsv%2Egz

Create experiment file (``experiments.csv``):

.. code-block:: text

    Condition,Replicate,DNA_BC_F,DNA_UMI,DNA_BC_R,RNA_BC_F,RNA_UMI,RNA_BC_R
    neurons,1,SRR37425536_1.fastq.gz,SRR37425536_2.fastq.gz,SRR37425536_3.fastq.gz,SRR37425544_1.fastq.gz,SRR37425544_2.fastq.gz,SRR37425544_3.fastq.gz
    neurons,2,SRR37425537_1.fastq.gz,SRR37425537_2.fastq.gz,SRR37425537_3.fastq.gz,SRR37425545_1.fastq.gz,SRR37425545_2.fastq.gz,SRR37425545_3.fastq.gz
    neurons,3,SRR37425531_1.fastq.gz,SRR37425531_2.fastq.gz,SRR37425531_3.fastq.gz,SRR37425546_1.fastq.gz,SRR37425546_2.fastq.gz,SRR37425546_3.fastq.gz
    neurons,4,SRR37425532_1.fastq.gz,SRR37425532_2.fastq.gz,SRR37425532_3.fastq.gz,SRR37425533_1.fastq.gz,SRR37425533_2.fastq.gz,SRR37425533_3.fastq.gz


Configure Experiment Workflow
-----------------------------

Create a config file (e.g. ``config_experiment.yaml``):

.. code-block:: yaml

    ---
    version: "0.7.0"
    experiments:
      GSE325256:
        split_number: 30 # increases speed due to large input data
        bc_length: 15
        umi_length: 16
        data_folder: data/experiment
        experiment_file: experiments.csv
        assignments:
          GSE325256BWAAssignment:
            type: file
            assignment_file: results/assignment/GSE325256BWA/assignment_barcodes.default.tsv.gz
          GSE325256CIGARAssignment:
            type: file
            assignment_file: results/assignment/GSE325256CIGAR/assignment_barcodes.default.tsv.gz
          GSE325256Assignment:
            type: file
            assignment_file: results/assignment/GSE325256/assignment_barcodes.default.tsv.gz
          GEOAssignment:
            type: file
            assignment_file: data/GSM9640449_VariantMPRAbarcodeDictionary.tsv.gz
        configs:
          default: {}


Run Experiment Workflow
-----------------------

.. code-block:: bash

    conda activate mprasnakeflow
    snakemake -c 1 --sdm conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile config_experiment.yaml -n -q --set-threads assignment_mapping_bbmap=10 --resources mem_mb=60000
    snakemake -c 30 --sdm conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile config_experiment.yaml -q --set-threads assignment_mapping_bbmap=10 --resources mem_mb=60000


Results
=======

Assignment results are written to :code:`results/assignment/GSE325256`, :code:`results/assignment/GSE325256BWA`, :code:`results/assignment/GSE325256CIGAR`, and :code:`results/assignment/GSE325256StrandSensitive`. The final assignment files and QC reports are placed in those directories. The library is very deeply sequenced and we get a median of over 1000 barcodes per oligo. Due to the challenging design (e.g. contains one base deletion but not within but at the end of the designed oligo) we also do not get the whole library assigned. Only around 70%. We can increase the assignment allowing a lower fraction percentage of reads assigned to the top hit (e.g. 0.51 instead of 0.75) but that also increases the risk of false positive assignments. Also the customized MPRAsnakeflow mapping with ``bwa`` disabling ``MAPQ`` with ``min_mapping_quality: 0`` but allwoing only specific CIGAR strings ``cigar_filter_regex: 208M|210M|212M|214M|216M|218M|219M`` works sightly better.

To reduce false positives, we also added the strand-sensitive mapping option (``GSE325256StrandSensitive`` using the default ``bbmap`` approach). We assign on average half of the barcodes (~500 per oligo). But it can be that we reduce the false positive rate.

The Experiment QC report later can help to make an informed decision on the best assignment parameters for this dataset.

You can have a look at the assignment QC reports here:

.. list-table:: Assignment QC reports for GSE325256
   :header-rows: 1

   * - Assignment mode
     - Fraction
     - QC report
   * - bbmap
     - 0.75
     - `Open report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE325256/GSE325256.assignment.qc_report.bbmap.default.html>`_
   * - bbmap
     - 0.51
     - `Open report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE325256/GSE325256.assignment.qc_report.bbmap.lowFraction.html>`_
   * - bwa-additional-filtering
     - 0.75
     - `Open report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE325256/GSE325256.assignment.qc_report.bwa_additional_filtering.default.html>`_
   * - bwa-additional-filtering
     - 0.51
     - `Open report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE325256/GSE325256.assignment.qc_report.bwa_additional_filtering.lowFraction.html>`_
   * - bwa with CIGAR filter
     - 0.75
     - `Open report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE325256/GSE325256.assignment.qc_report.bwa.default.html>`_
   * - bwa with CIGAR filter
     - 0.51
     - `Open report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE325256/GSE325256.assignment.qc_report.bwa.lowFraction.html>`_
   * - bbmap strand-sensitive mode
     - 0.75
     - `Open report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE325256/GSE325256.assignment.qc_report.bwa.default.html>`_
   * - bbmap strand-sensitive mode
     - 0.51
     - `Open report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE325256/GSE325256.assignment.qc_report.strand-sensitive.lowFraction.html>`_


Experiment results are written to :code:`results/experiments/GSE325256`. The count outputs use the published GEO assignment as well as the assignments generated in this workflow. We only used the mapping with the default fraction (0.75) for the experiment results.

The results show that in comparison to their assignment file we reached a similar pearosn correlation (around 0.67) for activity. Interestingly the strand-sensitive mode increased the correlsion to 0.72 suggesting that a high number of barcodes are false positives and the data might benefit from a strict mappinng strategy.


.. list-table:: Experiment QC reports for GSE325256
   :header-rows: 1

   * - Assignment mode
     - QC report
   * - GEO assignment
     - `Open report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE325256/GSE325256.experiment.qc_report.geo.html>`_
   * - bbmap
     - `Open report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE325256/GSE325256.experiment.qc_report.bbmap.default.html>`_
   * - bwa-additional-filtering
     - `Open report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE325256/GSE325256.experiment.qc_report.bwa_additional_filtering.default.html>`_
   * - bwa with CIGAR filter
     - `Open report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE325256/GSE325256.experiment.qc_report.bwa.default.html>`_
   * - bwa strand-sensitive mode
     - `Open report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE325256/GSE325256.experiment.qc_report.strand-sensitive.default.html>`_
