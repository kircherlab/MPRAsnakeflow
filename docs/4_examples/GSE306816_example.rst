.. _GSE306816 example:

.. role:: bash(code)
    :language: bash

===================================
GSE306816 (Zhang et al. dpSTR MPRA)
===================================

This example runs the assignment and experiment workflows on data from the preprint by `Zhang et al. Systematic Evaluation of the Impact of Promoter Proximal Short Tandem Repeats on Expression. bioRxiv (2025). <https://doi.org/10.1101/2025.09.14.676153>`_. The data were published in `GEO:GSE306816 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306816>`_.

The authors used a custom processing pipeline, available on `GitHub <https://github.com/gymreklab/str_mpra_design>`_, but here we run the analysis in MPRAsnakeflow. We only use their deep perturbation STR data (dpSTR). They also generated two additional MPRA datasets, h1STR and h2STR.

.. note:: This work measures perturbations of short tandem repeats, which are difficult to sequence because of homopolymer stretches and repeats. Zhang et al. used a specific stutter-correction method during assignment. That method is highly specific to their data and cannot be generalized. Even without stutter correction, we obtain assignment results comparable to those reported in their manuscript. The counting/experiment workflow is very similar, and we can also use their published assignment files.


Prerequisites
=============

This example depends on the following data and software:


Installation of MPRAsnakeflow
------------------------------

Please install conda, set up the MPRAsnakeflow environment, and clone the current MPRAsnakeflow master branch. You can find more help under :ref:`Installation`.

In addition, install `sra-toolkit` to download the data from GEO (for example via conda):

.. code-block:: bash

    conda install -c bioconda sra-tools

Design file
-----------

We need the design file and must modify it by trimming non-sequenced sequence parts and deduplicating the oligos used for assignment.

.. code-block:: bash

    mkdir -p data/assignment

    wget -O data/assignment/uber_seq.fa https://raw.githubusercontent.com/gymreklab/str-mpra/refs/heads/master/design-uber/uber_seq.fa

    # Use only the first 135 bp of each oligo and remove duplicates.
    # The full oligo design is 231 bp, but only 135 bp are sequenced.
    # Duplicates can exist because one oligo can map to multiple barcodes.
    awk -v N=135 '
    function flush(   s,key) {
        if (hdr == "") return
        gsub(/[ \t\r\n]/, "", seq)
        s = substr(seq, 1, N)
        if (s == "") return

        if (!(s in seen)) {
        seen[s] = 1
        order[++n] = s
        headers[s] = hdr
        } else {
        headers[s] = headers[s] "|" hdr
        }
        seq = ""; hdr = ""
    }

    /^>/ { flush(); hdr = substr($0,2); next }
    { seq = seq $0 }
    END { flush();
            for (i=1; i<=n; i++) {
            key = order[i]
            print ">" headers[key]
            print key
            }
        }
    ' data/assignment/uber_seq.fa > data/assignment/uber_seq_135pb_dedup.fa


Reads assignment data
----------------------

There are two sequencing sets: one generated on an Illumina sequencer and one on an Element sequencer. Here we process both together in parallel (SRA IDs: `SRR35184122` and `SRR35184121`).

.. code-block:: bash

    mkdir -p data/assignment
    cd data/assignment
    prefetch --max-size 130GB SRR35184122 SRR35184121
    fastq-dump --gzip --split-files SRR35184122
    fastq-dump --gzip --split-files SRR35184121

    cd ../../

.. note:: Please be sure that all files are downloaded completely without errors!

With

.. code-block:: bash

    tree data


The folder should look like this:

.. code-block:: text

    data
    └── assignment
       ├── SRR35184121_1.fastq.gz
       ├── SRR35184121_2.fastq.gz
       ├── SRR35184122_1.fastq.gz
       ├── SRR35184122_2.fastq.gz
       ├── uber_seq.fa
       └── uber_seq_135pb_dedup.fa

Reads count data
----------------

We have three replicates of DNA and RNA counts each. These data must be downloaded.

.. list-table::
    :header-rows: 1

    * - Replicate
        - gDNA
        - cDNA
    * - 1
        - SRR35184099
        - SRR35184102
    * - 2
        - SRR35184098
        - SRR35184101
    * - 3
        - SRR35184097
        - SRR35184100

.. code-block:: bash

    mkdir -p data/experiment
    cd data/experiment
    # download data
    prefetch --max-size 30GB SRR35184102 SRR35184101 SRR35184100 SRR35184099 SRR35184098 SRR35184097;
    fastq-dump --gzip --split-files SRR35184102 && \
    fastq-dump --gzip --split-files SRR35184101 && \
    fastq-dump --gzip --split-files SRR35184100 && \
    fastq-dump --gzip --split-files SRR35184099 && \
    fastq-dump --gzip --split-files SRR35184098 && \
    fastq-dump --gzip --split-files SRR35184097

    cd ../../

With

.. code-block:: bash

    tree data


The folder should look like this:

.. code-block:: text

    data
    ├── assignment
    │   ├── SRR35184121_1.fastq.gz
    │   ├── SRR35184121_2.fastq.gz
    │   ├── SRR35184122_1.fastq.gz
    │   ├── SRR35184122_2.fastq.gz
    │   ├── uber_seq.fa
    │   └── uber_seq_135pb_dedup.fa
    └── experiment
        ├── SRR35184097_1.fastq.gz
        ├── SRR35184098_1.fastq.gz
        ├── SRR35184099_1.fastq.gz
        ├── SRR35184100_1.fastq.gz
        ├── SRR35184101_1.fastq.gz
        └── SRR35184102_1.fastq.gz


Their assignment data
---------------------

From the GEO dataset, we can also download the assignment files they used for counting (Element and Illumina separately). We use these files for the experiment workflow. In addition, we run the assignment workflow on the raw FASTQ files to show that we obtain comparable results.

.. code-block:: bash

    mkdir -p data
    cd data

    wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9209nnn/GSM9209852/suppl/GSM9209852%5FdpSTR%5Fassociation%5FIllumina%5FBC%5FSTRs%2Etsv%2Egz
    wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9209nnn/GSM9209853/suppl/GSM9209853%5FdpSTR%5Fassociation%5FElement%5FBC%5FSTRs%2Etsv%2Egz

    zcat GSM9209852_dpSTR_association_Illumina_BC_STRs.tsv.gz  | tail -n+4 | gzip -c > GSM9209852_dpSTR_association_Illumina_BC_STRs_raw.tsv.gz
    zcat GSM9209853_dpSTR_association_Element_BC_STRs.tsv.gz | tail -n+4 | gzip -c > GSM9209853_dpSTR_association_Element_BC_STRs_raw.tsv.gz

    cd ../


MPRAsnakeflow
=============

We run assignment and count workflows together. It is also possible to run them separately using different config files. In that case, use assignment ``fromFile`` instead of ``fromConfig``.

First, define the config file and the experiment CSV file to map DNA/RNA counts to the correct replicates.

Important details for the config file: the oligo design is 135 bp and barcodes are 20 bp long. Remove the 5-bp sequence (``TCTAG``) at the 3-prime end of the barcode read.

For counting, reads are longer than the 20-bp barcode. There is also a 3-prime adapter sequence that we remove using cutadapt.

Create config files
-------------------

.. code-block:: bash

    cat << 'EOF' >  config.yaml
    ---
    version: "0.7"
    assignments:
        dpSTR:
            bc_length: 20
            adapters:
                BC:
                    three_prime:
                        - TCTAG
            alignment_tool:
                split_number: 30
                tool: bbmap
                configs:
                    min_mapping_quality: 30
            design_check:
                sequence_start: 1
                sequence_length: 135
            FWD:
                - data/assignment/SRR35184121_2.fastq.gz
                - data/assignment/SRR35184122_2.fastq.gz
            BC:
                - data/assignment/SRR35184121_1.fastq.gz
                - data/assignment/SRR35184122_1.fastq.gz
            design_file: data/assignment/uber_seq_135pb_dedup.fa
            configs:
                default: {}
    experiments:
        dpSTR:
            bc_length: 20
            adapters:
                FWD:
                    three_prime:
                        - TCTAGAGGTTCGTCGACGCGATTATTATCATTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAA
            data_folder: data/experiment
            experiment_file: experiments.csv
            assignments:
                dpSTRAssignment:
                    type: config
                    assignment_name: dpSTR
                    assignment_config: default
                GSM9209852Illumina:
                    type: file
                    assignment_file: data/GSM9209852_dpSTR_association_Illumina_BC_STRs_raw.tsv.gz
                GSM9209853Element:
                    type: file
                    assignment_file: data/GSM9209853_dpSTR_association_Element_BC_STRs_raw.tsv.gz
            configs:
                default: {}
    EOF

Create the :code:`experiments.csv` file to map DNA/RNA counts to replicates. The experiment file is a simple CSV file with the following content:

.. code-block:: bash

    cat << 'EOF' >  experiments.csv
    Condition,Replicate,DNA_BC_F,RNA_BC_F
    dpSTRHelaRNHWT,1,SRR35184099_1.fastq.gz,SRR35184102_1.fastq.gz
    dpSTRHelaRNHWT,2,SRR35184098_1.fastq.gz,SRR35184101_1.fastq.gz
    dpSTRHelaRNHWT,3,SRR35184097_1.fastq.gz,SRR35184100_1.fastq.gz
    EOF


Run snakemake
-------------

Now we are ready to run MPRAsnakeflow. We run this example on one node with 60 GB memory and 30 cores.

We run the pipeline directly in the working folder. The MPRAsnakeflow workflow can be located in a different directory. Here we assume :code:`/home/user/MPRAsnakeflow`.

First, do a dry run with snakemake using :code:`-n`:

.. code-block:: bash

    conda activate mprasnakeflow
    snakemake -c 1 --sdm apptainer conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile config.yaml -n --quiet rules

You should see a list of rules that will be executed. Example summary:

.. code-block:: text

    Job stats:
    job                                                                        count
    -----------------------------------------------------------------------  -------
    assignment_check_design                                                        1
    assignment_preprocessing_adapter_remove                                        1
    experiment_assigned_counts_filterAssignment                                    3
    experiment_preprocessing_trim_reads                                            6
    assignment_fastq_split                                                         2
    experiment_counts_onlyFWD_raw_counts                                           6
    assignment_attach_idx                                                         30
    experiment_counts_filter_counts                                                6
    experiment_statistic_counts_BC_in_RNA_DNA                                      6
    experiment_statistic_counts_table                                             12
    assignment_mapping_bbmap                                                      30
    experiment_counts_final_counts                                                 6
    experiment_statistic_counts_BC_in_RNA_DNA_merge                                2
    experiment_statistic_counts_frequent_umis                                      6
    experiment_statistic_counts_stats_merge                                        2
    assignment_collect                                                             1
    assignment_mapping_bbmap_getBCs                                               30
    experiment_assigned_counts_assignBarcodes                                     18
    experiment_counts_dna_rna_merge_counts                                        12
    experiment_statistic_bc_overlap_run                                            8
    experiment_statistic_counts_barcode_base_composition                           6
    experiment_statistic_counts_final                                              2
    assignment_collectBCs                                                          1
    assignment_idx_bam                                                             1
    experiment_assigned_counts_dna_rna_merge                                       9
    experiment_statistic_assigned_counts_combine_BC_assignment_stats_helper        3
    experiment_statistic_bc_overlap_combine_counts                                 1
    experiment_statistic_correlation_bc_counts                                     4
    experiment_statistic_correlation_bc_counts_hist                                4
    assignment_filter                                                              1
    assignment_flagstat                                                            1
    assignment_statistic_totalCounts                                               1
    experiment_assigned_counts_combine_replicates_barcode_output                   3
    experiment_assigned_counts_make_master_tables                                  3
    experiment_statistic_assigned_counts_combine_BC_assignment_stats               3
    experiment_statistic_assigned_counts_combine_stats_dna_rna_merge               3
    experiment_statistic_bc_overlap_combine_assigned_counts                        3
    experiment_statistic_correlation_calculate                                     3
    experiment_statistic_correlation_combine_bc_raw                                1
    experiment_statistic_correlation_hist_box_plots                                3
    assignment_statistic_assignedCounts                                            1
    assignment_statistic_assignment                                                1
    assignment_statistic_quality_metric                                            1
    experiment_assigned_counts_combine_replicates                                  6
    experiment_assigned_counts_copy_final_all_files                                3
    experiment_assigned_counts_copy_final_thresh_files                             3
    experiment_statistic_assigned_counts_combine_stats_dna_rna_merge_all           3
    experiment_statistic_correlation_combine_bc_assigned                           3
    experiment_statistic_correlation_combine_oligo                                 3
    experiment_statistic_quality_metric                                            3
    qc_report_assoc                                                                1
    qc_report_count                                                                3
    all                                                                            1
    total                                                                        276


If the dry run finishes without errors, run the workflow. We use a machine with 30 threads/cores and 60 GB memory. Therefore, :code:`split_number` is set to 30 for parallelization. We also use 10 threads for mapping (bbmap), and snakemake ensures that no more than 30 threads are used in total.

.. code-block:: bash

    snakemake -c 30 --sdm conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile config.yaml -q --set-threads assignment_mapping_bbmap=10 --resources mem_mb=60000


.. note:: Please adapt your setup when running in a cluster environment. An example SLURM profile is available in MPRAsnakeflow under :code:`profiles/default/config.yaml`. You can use it with snakemake via :code:`--workflow-profile $PIPELINE/profiles/default`, but adjust it first, especially the :code:`slurm_partition` setting.

Results
-------

For assignment, all output files are written to :code:`results/assignment/dpSTR`. The final assignment file is :code:`results/assignment/dpSTR/assignment_barcodes.default.tsv.gz`. You should also inspect the QC report: :code:`results/assignment/dpSTR/qc_report.default.html`.

An example assignment QC report is available here: `Example assignment QC report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE306816/GSE306816.assignment.qc_report.default.html>`_.


For experiments, all output files are written to :code:`results/experiments/dpSTR`.

Final count files:

* :code:`results/experiments/dpSTR/reporter_experiment.barcode.dpSTRHelaRNHWT.dpSTRAssignment.default.all.tsv.gz` (assignment generated in this workflow)
* :code:`results/experiments/dpSTR/reporter_experiment.barcode.dpSTRHelaRNHWT.GSM9209852Illumina.default.all.tsv.gz` (published Illumina assignment)
* :code:`results/experiments/dpSTR/reporter_experiment.barcode.dpSTRHelaRNHWT.GSM9209853Element.default.all.tsv.gz` (published Element assignment)

You should also inspect the QC reports, for example :code:`results/experiments/dpSTR/qc_report.dpSTRHelaRNHWT.GSM9209852Illumina.default.html`.

Example QC reports are available here: `Experiment QC report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE306816/GSE306816.qc_report.dpSTRHelaRNHWT.dpSTRAssignment.default.html>`_, `Illumina QC report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE306816/GSE306816.qc_report.dpSTRHelaRNHWT.GSM9209852Illumina.default.html>`_, and `Element QC report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE306816/GSE306816.qc_report.dpSTRHelaRNHWT.GSM9209853Element.default.html>`_.
