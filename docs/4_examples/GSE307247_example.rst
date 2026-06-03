.. _GSE307247 example:

.. role:: bash(code)
    :language: bash

=======================
GSE307247 Koplik et al.
=======================

This example runs the experiment workflow on data from `Koplik et al. Massively parallel assay of human splice variants reveals cis-regulatory drivers of disease-associated and cell type-specific splicing regulation. bioRxiv. (2025). <https://www.biorxiv.org/content/10.1101/2025.10.12.681955v1.full#sec-19>`_. The data were published in `GEO:GSE307247 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307247>`_.

The authors used a custom processing pipeline, available on `GitHub <https://github.com/skoplik/ESL_MPRA_2025>`_, but here we run part of the analysis in MPRAsnakeflow.

.. note:: This experiment represents a significant departure from a traditional MPRA experiment. Specifically, barcodes are used to associate reads with reporter constructs, but the ultimate quantification relies on counting exon skipping or inclusion events. The reference sequences consist of specific exons and introns with variable flanking sequences. The highly customized steps are mostly bash and Python scripts, so they can be readily integrated into a Snakemake pipeline. Here, we merely demonstrate how a user could intersect this pipeline with MPRAsnakeflow for count quantification without actually integrating the entire workflow.


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

We need the design file and must modify it by parsing the csv file, removing the header, and reformatting for use with MPRAsnakeflow.

.. code-block:: bash

    mkdir -p data/Koplik

    wget -O data/Koplik/Koplik.barcode.dictionary.fasta.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE307nnn/GSE307247/suppl/GSE307247%5FESL%5Fconcat%5F2023%5F09%5F19%5Fsubsampleparams%5Fd1c%5Fms75%5Fshorter3p%5Fiterate%5Fmincov5%5Freference%2Efasta%2Egz

    zcat data/Koplik/Koplik.barcode.dictionary.fasta.gz | \
        awk '/^>/ {if (N>0) printf "\n"; printf "%s\t", $0; N++; next} {printf "%s", $0} END {if (N>0) printf "\n"}' | \
        sed 's/^>//g' | \
        sed 's/ /_/g' | \
        awk '{print substr($2, length($2)-19, 20)"\t"$1}' | \
        cut -f1 | \
        rev | \
        tr ACGTN TGCAN \
        > data/Koplik/Koplik.bacodes.temp

    zcat data/Koplik/Koplik.barcode.dictionary.fasta.gz | \
        awk '/^>/ {if (N>0) printf "\n"; printf "%s\t", $0; N++; next} {printf "%s", $0} END {if (N>0) printf "\n"}' | \
        sed 's/^>//g' | \
        sed 's/ /_/g' | \
        awk '{print substr($2, length($2)-19, 20)"\t"$1}' | \
        cut -f2 \
        > data/Koplik/Koplik.names.temp

    paste data/Koplik/Koplik.barcodes.temp data/Koplik/Koplik.names.temp | \
        gzip \
        > data/Koplik/Koplik.barcode.dictionary.tsv.gz

    rm data/Koplik/Koplik.bacodes.temp data/Koplik/Koplik.names.temp

Read experiment data
--------------------

There is only one set of sequencing data for this experiment, and we are selecting all experiments (HEK293, HeLa, HMC3, K562, MCF7), each having two replicates. Each replicate has two fastq files associated with it, which are combined and compressed for convenience and compatibility with MPRAsnakeflow, which expects compressed fastq files. Raw sequencing reads are obtained from GEO and processed together in a list of accessions. Note the need for the --include-technical flag to ensure obtaining the necessary sequencing files containing UMIs for deduplication.

.. code-block:: bash

    for i in {192..202}; do echo SRR35247${i}; done > GSE307247_Acc_List.txt

    prefetch --option-file GSE307247_Acc_List.txt

    while read -r acc; do
        fasterq-dump ${acc} --include-technical --split-files
    done < GSE307247_Acc_List.txt

    gzip -c SRR35247192_2.fastq > data/Koplik/SRR35247192_2.fastq.gz
    gzip -c SRR35247193_3.fastq > data/Koplik/SRR35247193_3.fastq.gz
    gzip -c SRR35247193_4.fastq > data/Koplik/SRR35247193_4.fastq.gz
    gzip -c SRR35247194_3.fastq > data/Koplik/SRR35247194_3.fastq.gz
    gzip -c SRR35247194_4.fastq > data/Koplik/SRR35247194_4.fastq.gz
    gzip -c SRR35247195_3.fastq > data/Koplik/SRR35247195_3.fastq.gz
    gzip -c SRR35247195_4.fastq > data/Koplik/SRR35247195_4.fastq.gz
    gzip -c SRR35247196_3.fastq > data/Koplik/SRR35247196_3.fastq.gz
    gzip -c SRR35247196_4.fastq > data/Koplik/SRR35247196_4.fastq.gz
    gzip -c SRR35247197_3.fastq > data/Koplik/SRR35247197_3.fastq.gz
    gzip -c SRR35247197_4.fastq > data/Koplik/SRR35247197_4.fastq.gz
    gzip -c SRR35247198_3.fastq > data/Koplik/SRR35247198_3.fastq.gz
    gzip -c SRR35247198_4.fastq > data/Koplik/SRR35247198_4.fastq.gz
    gzip -c SRR35247199_3.fastq > data/Koplik/SRR35247199_3.fastq.gz
    gzip -c SRR35247199_4.fastq > data/Koplik/SRR35247199_4.fastq.gz
    gzip -c SRR35247200_3.fastq > data/Koplik/SRR35247200_3.fastq.gz
    gzip -c SRR35247200_4.fastq > data/Koplik/SRR35247200_4.fastq.gz
    gzip -c SRR35247201_3.fastq > data/Koplik/SRR35247201_3.fastq.gz
    gzip -c SRR35247201_4.fastq > data/Koplik/SRR35247201_4.fastq.gz
    gzip -c SRR35247202_3.fastq > data/Koplik/SRR35247202_3.fastq.gz
    gzip -c SRR35247202_4.fastq > data/Koplik/SRR35247202_4.fastq.gz

.. note:: Please be sure that all files are downloaded completely without errors!

With

.. code-block:: bash

    tree data


The folder should look like this:

.. code-block:: text

    data
    └── Koplik
       ├── Koplik.barcode.dictionary.fasta.gz
       ├── Koplik.barcode.dictionary.tsv.gz
        ├── SRR35247192_2.fastq.gz
        ├── SRR35247193_3.fastq.gz
        ├── SRR35247193_4.fastq.gz
        ├── SRR35247194_3.fastq.gz
        ├── SRR35247194_4.fastq.gz
        ├── SRR35247195_3.fastq.gz
        ├── SRR35247195_4.fastq.gz
        ├── SRR35247196_3.fastq.gz
        ├── SRR35247196_4.fastq.gz
        ├── SRR35247197_3.fastq.gz
        ├── SRR35247197_4.fastq.gz
        ├── SRR35247198_3.fastq.gz
        ├── SRR35247198_4.fastq.gz
        ├── SRR35247199_3.fastq.gz
        ├── SRR35247199_4.fastq.gz
        ├── SRR35247200_3.fastq.gz
        ├── SRR35247200_4.fastq.gz
        ├── SRR35247201_3.fastq.gz
        ├── SRR35247201_4.fastq.gz
        ├── SRR35247202_3.fastq.gz
        └── SRR35247202_4.fastq.gz

    1 directory, 23 files


Their assignment data
---------------------

From the GEO dataset, we can also download the counts reported by the authors for comparison to what we obtain using MPRAsnakeflow and show that we obtain comparable results.

.. code-block:: bash

    wget -O data/Koplik/Koplik.counts.csv.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE307nnn/GSE307247/suppl/GSE307247%5FProcessed%5FPSIs%5FAll%5FCells%2Ecsv%2Egz


MPRAsnakeflow
=============

We run only the count workflow. Note, we use assignment ``fromFile`` instead of ``fromConfig`` using the assignment file we generated from the supplied barcode map.

First, define the config file and the experiment CSV file to map DNA/RNA counts to the correct replicates.

Important details for the config file: the barcodes are 20 bp long. Deduplication using UMIs is performed, so it is critical to supply umi_length, 8 bp in this case. In addition, note that the default config settings are not used here. Due to the nature of the experiment, the barcode threshold needs to be set to 1 to avoid filtering out almost all barcodes.

Create config files
-------------------

.. code-block:: bash

    cat << 'EOF' >  config/Koplik.yaml
    ---
    version: "0.6.5"
    experiments:
        Koplik:
            bc_length: 20
            umi_length: 8
            data_folder: data/Koplik
            experiment_file: resources/Koplik.csv
            demultiplex: true
            assignments:
                fromFile:
                    type: file
                    assignment_file: data/Koplik/Koplik.barcode.dictionary.tsv.gz
            configs:
                bcone:
                    filter:
                        bc_threshold: 1
    EOF

Create the :code:`experiments.csv` file to map DNA/RNA counts to replicates. The experiment file is a simple CSV file with the following content:

.. code-block:: bash

    cat << 'EOF' >  resources/Koplik.csv
    Condition,Replicate,DNA_BC_F,RNA_BC_F,RNA_UMI
    HEK293,1,SRR35247192_2.fastq.gz,SRR35247202_4.fastq.gz,SRR35247202_3.fastq.gz
    HEK293,2,SRR35247192_2.fastq.gz,SRR35247201_4.fastq.gz,SRR35247201_3.fastq.gz
    HeLa,1,SRR35247192_2.fastq.gz,SRR35247200_4.fastq.gz,SRR35247200_3.fastq.gz
    HeLa,2,SRR35247192_2.fastq.gz,SRR35247199_4.fastq.gz,SRR35247199_3.fastq.gz
    HMC3,1,SRR35247192_2.fastq.gz,SRR35247198_4.fastq.gz,SRR35247198_3.fastq.gz
    HMC3,2,SRR35247192_2.fastq.gz,SRR35247197_4.fastq.gz,SRR35247197_3.fastq.gz
    K562,1,SRR35247192_2.fastq.gz,SRR35247196_4.fastq.gz,SRR35247196_3.fastq.gz
    K562,2,SRR35247192_2.fastq.gz,SRR35247195_4.fastq.gz,SRR35247195_3.fastq.gz
    MCF7,1,SRR35247192_2.fastq.gz,SRR35247194_4.fastq.gz,SRR35247194_3.fastq.gz
    MCF7,2,SRR35247192_2.fastq.gz,SRR35247193_4.fastq.gz,SRR35247193_3.fastq.gz
    EOF


Run snakemake
-------------

Now we are ready to run MPRAsnakeflow. We run this example on a HPC cluster using the slurm executor plugin for Snakemake

We run the pipeline directly in the working folder. The MPRAsnakeflow workflow can be located in a different directory. Here we assume :code:`/home/user/MPRAsnakeflow`.

First, do a dry run with snakemake using :code:`-n`:

.. code-block:: bash

    conda activate mprasnakeflow
    snakemake --software-deployment-method conda --executor slurm --jobs 24 --configfile config/Koplik.yaml --workflow-profile profiles/default -n --quiet rules

You should see a list of rules that will be executed. Example summary:

.. code-block:: text

    Job stats:
    job                                                                        count
    -----------------------------------------------------------------------  -------
    all                                                                            1
    experiment_assigned_counts_assignBarcodes                                     20
    experiment_assigned_counts_combine_replicates                                 10
    experiment_assigned_counts_combine_replicates_barcode_output                   5
    experiment_assigned_counts_copy_final_all_files                                5
    experiment_assigned_counts_copy_final_thresh_files                             5
    experiment_assigned_counts_dna_rna_merge                                      10
    experiment_assigned_counts_make_master_tables                                  5
    experiment_counts_dna_rna_merge_counts                                        20
    experiment_counts_filter_counts                                               15
    experiment_counts_final_counts                                                15
    experiment_counts_onlyFWDUMI_raw_counts                                       10
    experiment_counts_onlyFWD_raw_counts                                           5
    experiment_statistic_assigned_counts_combine_BC_assignment_stats               1
    experiment_statistic_assigned_counts_combine_BC_assignment_stats_helper        5
    experiment_statistic_assigned_counts_combine_stats_dna_rna_merge               5
    experiment_statistic_assigned_counts_combine_stats_dna_rna_merge_all           1
    experiment_statistic_bc_overlap_combine_assigned_counts                        1
    experiment_statistic_bc_overlap_combine_counts                                 1
    experiment_statistic_bc_overlap_run                                           20
    experiment_statistic_correlation_bc_counts                                    10
    experiment_statistic_correlation_bc_counts_hist                               10
    experiment_statistic_correlation_calculate                                     5
    experiment_statistic_correlation_combine_bc_assigned                           1
    experiment_statistic_correlation_combine_bc_raw                                1
    experiment_statistic_correlation_combine_oligo                                 1
    experiment_statistic_correlation_hist_box_plots                                5
    experiment_statistic_counts_BC_in_RNA_DNA                                     20
    experiment_statistic_counts_BC_in_RNA_DNA_merge                                2
    experiment_statistic_counts_barcode_base_composition                          15
    experiment_statistic_counts_final                                              2
    experiment_statistic_counts_frequent_umis                                     15
    experiment_statistic_counts_stats_merge                                        2
    experiment_statistic_counts_table                                             30
    experiment_statistic_quality_metric                                            5
    qc_report_count                                                                5
    total                                                                        289

If the dry run finishes without errors, run the workflow. We use a machine with 24 threads/cores. Therefore, :code:`split_number` is set to 24 for parallelization, and snakemake ensures that no more than 24 threads are used in total.

.. code-block:: bash

    sbatch -t 2- --mem 32G -o Koplik_count.out --wrap "snakemake --software-deployment-method conda --executor slurm --jobs 24 --configfile config/Koplik.yaml --workflow-profile profiles/default"


.. note:: An example SLURM profile is available in MPRAsnakeflow under :code:`profiles/default/config.yaml`. You can use it with snakemake via :code:`--workflow-profile $PIPELINE/profiles/default`, but adjust it first, especially the :code:`slurm_partition` setting.

Results
-------

For experiments, all output files are written to :code:`results/experiments/Koplik`.

Final count files:

* :code:`results/experiments/Koplik/reporter_experiment.barcode.HEK293.fromFile.bcone.all.tsv.gz` (counts generated in this workflow)
* :code:`results/experiments/Koplik/reporter_experiment.barcode.HeLa.fromFile.bcone.all.tsv.gz` (counts generated in this workflow)
* :code:`results/experiments/Koplik/reporter_experiment.barcode.HMC3.fromFile.bcone.all.tsv.gz` (counts generated in this workflow)
* :code:`results/experiments/Koplik/reporter_experiment.barcode.K562.fromFile.bcone.all.tsv.gz` (counts generated in this workflow)
* :code:`results/experiments/Koplik/reporter_experiment.barcode.MCF7.fromFile.bcone.all.tsv.gz` (counts generated in this workflow)

You should also inspect the QC reports, for example :code:`results/experiments/Koplik/qc_report.HEK293.fromFile.bcone.html`.
