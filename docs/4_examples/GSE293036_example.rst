.. _GSE293036 example:

.. role:: bash(code)
    :language: bash

=========
GSE293036
=========

This example runs the assignment and experiment workflows on data from `Granitto et al. Genome-wide discovery of multiple sclerosis genetic risk variant allelic regulatory activity. G3 Genes|Genomes|Genetics (2025) <https://doi.org/10.1093/g3journal/jkaf192>`_. The data were published in `GEO:GSE293036 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293036>`_.

The authors used a massively parallel reporter assay (MPRA) to identify multiple sclerosis (MS) genetic risk variants with allele-specific regulatory activity, testing oligos in three conditions: GM12878, MS-1, and MS-2.

.. note:: The authors shared their experimental design and data through GEO. For the barcode library, they aligned reads to a custom reference including both forward and reverse complement versions of the oligo sequences. For simplicity, we just included the forward versions in the design file, achieving similar oligo coverage (number of barcodes per oligo).

Prerequisites
=============

This example depends on the following data and software:


Installation of MPRAsnakeflow
------------------------------

Please install conda, set up the MPRAsnakeflow environment, and clone the current MPRAsnakeflow master branch. You can find more help under :ref:`Installation`.

In addition, install `sra-toolkit` and `openpyxl` to download and process the data (for example via conda):

.. code-block:: bash

    conda install -c bioconda sra-tools
    conda install -c conda-forge openpyxl

Design file
-----------

The oligo library sequences are provided in Supplementary Table 3 of the paper (``Supplementary_Table_3_G3-2025-406100.xlsx``, columns ``seq_id`` and ``twist_seq``). Download the supplementary data zip, extract the table, and convert it to a gzipped FASTA, filtering out the reverse-complement (``_RC``) entries:

.. code-block:: bash

    mkdir -p data
    cd data

    wget -O jkaf192_supplementary_data.zip \
        https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/g3journal/15/11/10.1093_g3journal_jkaf192/1/jkaf192_supplementary_data.zip
    unzip -j jkaf192_supplementary_data.zip Supplementary_Table_3_G3-2025-406100.xlsx

    python3 - << 'PYEOF'
    import openpyxl, gzip
    wb = openpyxl.load_workbook("Supplementary_Table_3_G3-2025-406100.xlsx", read_only=True)
    ws = wb.active
    with gzip.open("design.fa.gz", "wt") as fh:
        for row in ws.iter_rows(min_row=2, values_only=True):
            seq_id, twist_seq = row[0], row[1]
            if seq_id and not str(seq_id).endswith("_RC"):
                fh.write(f">{seq_id}\n{twist_seq}\n")
    PYEOF

    cd ..


Reads assignment data
----------------------

The barcode library sequencing (`GEO:GSM8874361 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8874361>`_) is available under SRA accession ``SRR32870551``. This paired-end run has the barcode in read 1 and the oligo sequence in read 2.

.. code-block:: bash

    mkdir -p data
    cd data
    prefetch SRR32870551
    fastq-dump --gzip --split-files SRR32870551 && \
    mv SRR32870551_1.fastq.gz SRR32870551_Barcode_1.fastq.gz && \
    mv SRR32870551_2.fastq.gz SRR32870551_Barcode_2.fastq.gz
    cd ..

.. note:: Please be sure that all files are downloaded completely without errors!

With

.. code-block:: bash

    tree data


The folder should look like this:

.. code-block:: text

    data
    ├── SRR32870551_Barcode_1.fastq.gz
    ├── SRR32870551_Barcode_2.fastq.gz
    └── design.fa.gz

Reads count data
----------------

In this example we run just one sequencing lane for two replicates of the GM12878 condition, using the corresponding control replicates as input.

.. list-table::
    :header-rows: 1

    * - Condition
      - Replicate
      - GSM
      - SRR
    * - Control
      - 2
      - GSM8874363
      - SRR32870549
    * - Control
      - 3
      - GSM8874364
      - SRR32870548
    * - GM12878
      - 2
      - GSM8874366
      - SRR32870540
    * - GM12878
      - 3
      - GSM8874367
      - SRR32870536

.. note:: Each GSM has multiple sequencing runs in SRA. This example uses only the first run (r1) of each replicate. For a complete analysis, download and include all runs per replicate.

.. code-block:: bash

    mkdir -p data
    cd data
    prefetch SRR32870549 SRR32870548 SRR32870540 SRR32870536
    fastq-dump --gzip --split-files SRR32870549 && mv SRR32870549_1.fastq.gz SRR32870549_Control.fastq.gz && \
    fastq-dump --gzip --split-files SRR32870548 && mv SRR32870548_1.fastq.gz SRR32870548_Control.fastq.gz && \
    fastq-dump --gzip --split-files SRR32870540 && mv SRR32870540_1.fastq.gz SRR32870540_GM12878.fastq.gz && \
    fastq-dump --gzip --split-files SRR32870536 && mv SRR32870536_1.fastq.gz SRR32870536_GM12878.fastq.gz
    cd ..

With

.. code-block:: bash

    tree data


The folder should look like this:

.. code-block:: text

    data
    ├── SRR32870536_GM12878.fastq.gz
    ├── SRR32870540_GM12878.fastq.gz
    ├── SRR32870548_Control.fastq.gz
    ├── SRR32870549_Control.fastq.gz
    ├── SRR32870551_Barcode_1.fastq.gz
    ├── SRR32870551_Barcode_2.fastq.gz
    └── design.fa.gz


MPRAsnakeflow
=============

We run assignment and count workflows together. It is also possible to run them separately using different config files. In that case, use assignment ``fromFile`` instead of ``fromConfig``.

First, define the config file and the experiment CSV file to map DNA/RNA counts to the correct replicates.

Important details for the config file: oligos are 200 bp long and barcodes are 20 bp. The linker sequence separating the barcode (read 1) from the oligo is :code:`TCTAGAGGTTCGTCGACGCGATCGCAGGAGCCGCAGTG`.

Create config files
-------------------

.. code-block:: bash

    cat << 'EOF' >  config.yaml
    ---
    version: "0.6.5"
    assignments:
        GSE293036Assignment:
            bc_length: 20
            BC_rev_comp: false
            linker: TCTAGAGGTTCGTCGACGCGATCGCAGGAGCCGCAGTG
            alignment_tool:
                split_number: 16
                tool: bbmap
                configs:
                    sequence_length: 200
                    alignment_start: 1
            FWD:
            - data/SRR32870551_Barcode_1.fastq.gz
            REV:
            - data/SRR32870551_Barcode_2.fastq.gz
            design_file: data/design.fa.gz
            configs:
                default: {}
    experiments:
        GSE293036Experiment:
            bc_length: 20
            split_number: 16
            data_folder: data
            experiment_file: experiments.csv
            assignments:
                GSE293036Assignment:
                    type: config
                    assignment_name: GSE293036Assignment
                    assignment_config: default
            configs:
                default:
                    filter:
                        bc_threshold: 30
    EOF

Create the :code:`experiments.csv` file to map DNA/RNA counts to replicates. The control samples serve as gDNA input (DNA_BC_F) for all conditions.

.. code-block:: bash

    cat << 'EOF' >  experiments.csv
    Condition,Replicate,DNA_BC_F,RNA_BC_F
    GM12878,3,SRR32870548_Control.fastq.gz,SRR32870536_GM12878.fastq.gz
    GM12878,2,SRR32870549_Control.fastq.gz,SRR32870540_GM12878.fastq.gz
    EOF


Run snakemake
-------------

Now we are ready to run MPRAsnakeflow. We run this example on one node with 60 GB memory and 16 cores.

We run the pipeline directly in the working folder. The MPRAsnakeflow workflow can be located in a different directory. Here we assume :code:`/work/${USER}/MPRAsnakeflow`.

First, do a dry run with snakemake using :code:`-n`:

.. code-block:: bash

    conda activate mprasnakeflow
    snakemake -c 1 --sdm apptainer conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile config.yaml -n --quiet rules

If the dry run finishes without errors, run the workflow. We use a machine with 16 threads/cores and 60 GB memory. Therefore, :code:`split_number` is set to 16 for parallelization. We also use 8 threads for mapping (bbmap), and snakemake ensures that no more than 16 threads are used in total.

.. code-block:: bash

    snakemake -c 16 --sdm apptainer conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --configfile config.yaml -q --set-threads assignment_mapping_bbmap=8 --resources mem_mb=60000


.. note:: Please adapt your setup when running in a cluster environment. An example SLURM profile is available in MPRAsnakeflow under :code:`profiles/default/config.yaml`. You can use it with snakemake via :code:`--workflow-profile $PIPELINE/profiles/default`, but adjust it first, especially the :code:`slurm_partition` setting.

Results
-------

For assignment, all output files are written to :code:`results/assignment/GSE293036Assignment`. The final assignment file is :code:`results/assignment/GSE293036Assignment/assignment_barcodes.default.tsv.gz`. You should also inspect the QC report: :code:`results/assignment/GSE293036Assignment/qc_report.default_GSE293036.html`.

An example assignment QC report is available here: `Example assignment QC report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE293036.assignment.qc_report.default.html>`_.


For experiments, all output files are written to :code:`results/experiments/GSE293036Experiment`.

Final count file:

* :code:`results/experiments/GSE293036Experiment/reporter_experiment.barcode.GM12878.GSE293036Assignment.default.all.tsv.gz`

You should also inspect the QC report: :code:`results/experiments/GSE293036Experiment/qc_report.GM12878.GSE293036Assignment.default.html`.

An example QC report is available here: `GM12878 QC report <https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAsnakeflow/blob/master/docs/4_examples/GSE293036.experiment.GM12878.qc_report.default.html>`_.
