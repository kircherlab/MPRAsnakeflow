.. _GSE271608 example:

.. role:: bash(code)
    :language: bash

=====================
GSE271608 Zahm et al.
=====================

This example runs the experiment workflow on data from `Zahm et al. A massively parallel reporter assay library to screen short synthetic promoters in mammalian cells. Nat Commun. (2024). <https://doi.org/10.1038/s41467-024-54502-9>`_. The data were published in `GEO:GSE271608 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271608>`_.

The authors used a custom processing pipeline, available on `GitHub <https://github.com/JGEnglishLab/TRE-MPRA-Pipeline>`_, but here we run part of the analysis in MPRAsnakeflow.

.. note:: The authors did not provide the fastq files to build the barcode dictionary so we demonstrate integrating the experiment workflow into an analysis where the equivalent of the association workflow is run separately.


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

    mkdir -p data/Zahm

    wget -O data/Zahm/Zahm.bacode.dictionary.csv.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271608/suppl/GSE271608%5FfinalBarcodeMap%2Ecsv%2Egz

    # Convert supplied csv file to format compatible with MPRAsnakeflow.
    zcat data/Zahm/Zahm.bacode.dictionary.csv.gz | \
        sed 's/,/\t/g' | \
        sed '1d' | \
        awk '{print $7"\t"$2"_"$3"_"$4"_"$5}' | \
        gzip > data/Zahm/Zahm.bacode.dictionary.tsv.gz


Read experiment data
--------------------

There is only one set of sequencing data for this experiment, and we are selecting two experiments (sf_19664 and sf_19919) that both have four replicates. Each replicate has two fastq files associated with it, which are combined and compressed for convenience and compatibility with MPRAsnakeflow, which expects compressed fastq files. Raw sequencing reads are obtained from GEO and processed together in a list of accessions.

.. code-block:: bash

    cd data/Zahm

    for i in {862..873} {876..879} {918..921}; do echo SRR29718${i}; done > GSE271608_Acc_List.txt

    prefetch --option-file GSE271608_Acc_List.txt

    while read acc; do
        fasterq-dump ${acc} --split-files
    done < GSE271608_Acc_List.txt

    cat SRR29718876.fastq SRR29718877.fastq | gzip > sf_19664.RNA1.fastq.gz
    cat SRR29718872.fastq SRR29718873.fastq | gzip > sf_19664.RNA2.fastq.gz
    cat SRR29718868.fastq SRR29718869.fastq | gzip > sf_19664.RNA3.fastq.gz
    cat SRR29718878.fastq SRR29718879.fastq | gzip > sf_19664.RNA4.fastq.gz
    cat SRR29718870_1.fastq SRR29718871_1.fastq | gzip > sf_19919.RNA1.fastq.gz
    cat SRR29718866_1.fastq SRR29718867_1.fastq | gzip > sf_19919.RNA2.fastq.gz
    cat SRR29718864_1.fastq SRR29718865_1.fastq | gzip > sf_19919.RNA3.fastq.gz
    cat SRR29718862_1.fastq SRR29718863_1.fastq | gzip > sf_19919.RNA4.fastq.gz
    cat SRR29718918_1.fastq SRR29718919_1.fastq | gzip > sf_19664.DNA.fastq.gz
    cat SRR29718920_1.fastq SRR29718921_1.fastq | gzip > sf_19919.DNA.fastq.gz

    cd ../../

.. note:: Please be sure that all files are downloaded completely without errors!

With

.. code-block:: bash

    tree data


The folder should look like this:

.. code-block:: text

    data
    └── Zahm
        ├── sf_19664.DNA.fastq.gz
        ├── sf_19664.RNA1.fastq.gz
        ├── sf_19664.RNA2.fastq.gz
        ├── sf_19664.RNA3.fastq.gz
        ├── sf_19664.RNA4.fastq.gz
        ├── sf_19919.DNA.fastq.gz
        ├── sf_19919.RNA1.fastq.gz
        ├── sf_19919.RNA2.fastq.gz
        ├── sf_19919.RNA3.fastq.gz
        ├── sf_19919.RNA4.fastq.gz
        ├── Zahm.bacode.dictionary.csv.gz
        └── Zahm.bacode.dictionary.tsv.gz

    1 directory, 12 files


Their assignment data
---------------------

From the GEO dataset, we can also download the counts reported by the authors for comparison to what we obtain using MPRAsnakeflow and show that we obtain comparable results.

.. code-block:: bash

    wget -O data/Zahm/Zahm.counts.csv.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271608/suppl/GSE271608%5Frpms%2Ecsv%2Egz


MPRAsnakeflow
=============

We run only the count workflow. Note, we use assignment ``fromFile`` instead of ``fromConfig`` using the assignment file we generated from the supplied barcode map.

First, define the config file and the experiment CSV file to map DNA/RNA counts to the correct replicates.

Important details for the config file: the barcodes are 24 bp long. Deduplication using UMIs is not done so umi_length is not used.

Create config files
-------------------

.. code-block:: bash

    cat << 'EOF' >  config/Zahm.yaml
    ---
    version: "0.6.5"
    experiments:
        Zahm:
            bc_length: 24
            umi_length: 16  # not used
            data_folder: data/Zahm
            experiment_file: resources/Zahm.csv
            demultiplex: false
            assignments:
                fromFile:
                    type: file
                    assignment_file: data/Zahm/Zahm.bacode.dictionary.tsv.gz
            configs:
                default: {} # name of an example filtering config
    EOF

Create the :code:`experiments.csv` file to map DNA/RNA counts to replicates. The experiment file is a simple CSV file with the following content:

.. code-block:: bash

    cat << 'EOF' >  resources/Zahm.csv
    Condition,Replicate,DNA_BC_F,RNA_BC_F
    sf19664,1,sf_19664.DNA.fastq.gz,sf_19664.RNA1.fastq.gz
    sf19664,2,sf_19664.DNA.fastq.gz,sf_19664.RNA2.fastq.gz
    sf19664,3,sf_19664.DNA.fastq.gz,sf_19664.RNA3.fastq.gz
    sf19664,4,sf_19664.DNA.fastq.gz,sf_19664.RNA4.fastq.gz
    sf19919,1,sf_19919.DNA.fastq.gz,sf_19919.RNA1.fastq.gz
    sf19919,2,sf_19919.DNA.fastq.gz,sf_19919.RNA2.fastq.gz
    sf19919,3,sf_19919.DNA.fastq.gz,sf_19919.RNA3.fastq.gz
    sf19919,4,sf_19919.DNA.fastq.gz,sf_19919.RNA4.fastq.gz
    EOF


Run snakemake
-------------

Now we are ready to run MPRAsnakeflow. We run this example on a HPC cluster using the slurm executor plugin for Snakemake

We run the pipeline directly in the working folder. The MPRAsnakeflow workflow can be located in a different directory. Here we assume :code:`/home/user/MPRAsnakeflow`.

First, do a dry run with snakemake using :code:`-n`:

.. code-block:: bash

    conda activate mprasnakeflow
    snakemake --software-deployment-method conda --executor slurm --jobs 24 --configfile config/Zahm.yaml --workflow-profile profiles/default -n --quiet rules

You should see a list of rules that will be executed. Example summary:

.. code-block:: text

    Job stats:
    job                                                                        count
    -----------------------------------------------------------------------  -------
    all                                                                            1
    experiment_assigned_counts_assignBarcodes                                     16
    experiment_assigned_counts_combine_replicates                                  4
    experiment_assigned_counts_combine_replicates_barcode_output                   2
    experiment_assigned_counts_copy_final_all_files                                2
    experiment_assigned_counts_copy_final_thresh_files                             2
    experiment_assigned_counts_dna_rna_merge                                       8
    experiment_assigned_counts_filterAssignment                                    1
    experiment_assigned_counts_make_master_tables                                  2
    experiment_counts_dna_rna_merge_counts                                        16
    experiment_counts_filter_counts                                               10
    experiment_counts_final_counts                                                10
    experiment_counts_onlyFWD_raw_counts                                          10
    experiment_statistic_assigned_counts_combine_BC_assignment_stats               1
    experiment_statistic_assigned_counts_combine_BC_assignment_stats_helper        2
    experiment_statistic_assigned_counts_combine_stats_dna_rna_merge               2
    experiment_statistic_assigned_counts_combine_stats_dna_rna_merge_all           1
    experiment_statistic_bc_overlap_combine_assigned_counts                        1
    experiment_statistic_bc_overlap_combine_counts                                 1
    experiment_statistic_bc_overlap_run                                            8
    experiment_statistic_correlation_bc_counts                                     4
    experiment_statistic_correlation_bc_counts_hist                                4
    experiment_statistic_correlation_calculate                                     2
    experiment_statistic_correlation_combine_bc_assigned                           1
    experiment_statistic_correlation_combine_bc_raw                                1
    experiment_statistic_correlation_combine_oligo                                 1
    experiment_statistic_correlation_hist_box_plots                                2
    experiment_statistic_counts_BC_in_RNA_DNA                                     16
    experiment_statistic_counts_BC_in_RNA_DNA_merge                                2
    experiment_statistic_counts_barcode_base_composition                          10
    experiment_statistic_counts_final                                              2
    experiment_statistic_counts_frequent_umis                                     10
    experiment_statistic_counts_stats_merge                                        2
    experiment_statistic_counts_table                                             20
    experiment_statistic_quality_metric                                            2
    qc_report_count                                                                2
    total                                                                        181

If the dry run finishes without errors, run the workflow. We use a machine with 24 threads/cores. Therefore, :code:`split_number` is set to 24 for parallelization, and snakemake ensures that no more than 24 threads are used in total.

.. code-block:: bash

    sbatch -t 2- --mem 32G -o Zahm_count.out --wrap "snakemake --software-deployment-method conda --executor slurm --jobs 24 --configfile config/Zahm.yaml --workflow-profile profiles/default"


.. note:: An example SLURM profile is available in MPRAsnakeflow under :code:`profiles/default/config.yaml`. You can use it with snakemake via :code:`--workflow-profile $PIPELINE/profiles/default`, but adjust it first, especially the :code:`slurm_partition` setting.

Results
-------

For experiments, all output files are written to :code:`results/experiments/Zahm`.

Final count files:

* :code:`results/experiments/Zahm/reporter_experiment.barcode.sf19664.fromFile.default.all.tsv.gz` (counts generated in this workflow)
* :code:`results/experiments/Zahm/reporter_experiment.barcode.sf19919.fromFile.default.all.tsv.gz` (counts generated in this workflow)

You should also inspect the QC reports, for example :code:`results/experiments/Zahm/qc_report.sf19664.fromFile.default.html`.
