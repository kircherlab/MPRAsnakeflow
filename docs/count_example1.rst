.. _Basic count workflow:

.. role:: bash(code)
      :language: bash

=====================
Basic Count workflow
=====================

This example runs the count workflow on 5'/5' WT MPRA data in the HEPG2 cell line from `Klein J., Agarwal, V., Keith, A., et al. 2019 <https://www.biorxiv.org/content/10.1101/576405v1.full.pdf>`_.

Prerequirements
======================

This example depends on the following data and software:


Installing MPRAsnakeflow
------------------------

Please install conda, the MPRAsnakeflow environment, and clone the actual ``MPRAsnakeflow`` master branch. You will find more help under :ref:`Installation`.

Producing an association (.tsv) file 
------------------------------------
This workflow requires a python dictionary of candidate regulatory sequence (CRS) mapped to their barcodes in a tab separated (.tsv) format. For this example the file can be generated using :ref:`Association example` or it can be found in `sample` folder in `MPRAsnakelfow <https://github.com/kircherlab/MPRAsnakeflow/>`_.

Alternatively, if the association file is in pickle (.pickle) format, you can convert the same file to .tsv format with the in-built function in MPRsnakeflow with the following code:

.. code-block:: bash
    
    conda activate mprasnakeflow
    python assignment_pickle_to_tsv.py --input <assignment_file in pickle format>


Design (.fa) file
-----------------

    File can be generated using the :ref:`Association example` or downloaded from the `sample` folder in `MPRAsnakelfow <https://github.com/kircherlab/MPRAsnakeflow/>`_.



Reads
-----

There is one condition (HEPG2) with three technical replicates. Each replicate contains a forward (barcode-forward), reverse (barcode-reverse), and index (unique molecular identifier) read for DNA and RNA. These data must be downloaded. All data is publically available on the short read archive (SRA). We will use SRA-toolkit to obtain the data.

.. note:: You need 9 GB disk space to download the data and upwards of 50 GB to proccess it!

.. code-block:: bash

    conda install sra-tools
    mkdir -p Count_Basic/data
    cd Count_Basic/data
    fastq-dump --gzip --split-files SRR10800881 SRR10800882 SRR10800883 SRR10800884 SRR10800885 SRR10800886
    cd ..

For large files and unstable internet connection we reccommend the comand ``prefetch`` from SRA tools before running ``fastq-dump``. This command is much smarter in warnings when something went wrong.

.. code-block:: bash

    conda install sra-tools
    cd Count_Basic/data
    prefetch SRR10800881 SRR10800882 SRR10800883 SRR10800884 SRR10800885 SRR10800886
    fastq-dump --gzip --split-files SRR10800986
    cd ..


.. note:: Please be sure that all files are downloaded completely without errors! Depending on your internet connection this can take a while. If you just want some data to run MPRAsnakeflow, you can just limit yourself to one condition and/or just one replicate.

The data folder view can be seen with the following command:

.. code-block:: bash

    tree data


The folder should look like this:

.. code-block:: text

    data

Here is an overview of the files:

.. csv-table:: HEPG2 data
   :header: "Condition", "GEO Accession", "SRA Accession", SRA Runs
   :widths: 40, 10, 10, 20

   "HEPG2-DNA-1: HEPG2 DNA replicate 1", GSM4237863, SRX7474781, "SRR10800881"
   "HEPG2-RNA-1: HEPG2 RNA replicate 1", GSM4237864, SRX7474782, "SRR10800882"
   "HEPG2-DNA-2: HEPG2 DNA replicate 2", GSM4237865, SRX7474783, "SRR10800883"
   "HEPG2-RNA-2: HEPG2 RNA replicate 2", GSM4237866, SRX7474784, "SRR10800884"
   "HEPG2-DNA-3: HEPG2 DNA replicate 3", GSM4237867, SRX7474785, "SRR10800885"
   "HEPG2-RNA-3: HEPG2 RNA replicate 3", GSM4237868, SRX7474786, "SRR10800886"



Run MPRAsnakeflow
=================

Now we are close to starting MPRAsnakeflow and count the number of barcodes. But before we need to generate an environment (.csv) file to tell snakemake the conditions, replicates and the corresponding reads.

Creating experiment.csv
---------------------------

Our experiment file looks exactly like this:

.. code-block:: text

    Condition,Replicate,DNA_BC_F,DNA_UMI,DNA_BC_R,RNA_BC_F,RNA_UMI,RNA_BC_R
    HEPG2,1,SRR10800881_1.fastq.gz,SRR10800881_2.fastq.gz,SRR10800881_3.fastq.gz,SRR10800882_1.fastq.gz,SRR10800882_2.fastq.gz,SRR10800882_3.fastq.gz
    HEPG2,2,SRR10800883_1.fastq.gz,SRR10800883_2.fastq.gz,SRR10800883_3.fastq.gz,SRR10800884_1.fastq.gz,SRR10800884_2.fastq.gz,SRR10800884_3.fastq.gz
    HEPG2,3,SRR10800885_1.fastq.gz,SRR10800885_2.fastq.gz,SRR10800885_3.fastq.gz,SRR10800886_1.fastq.gz,SRR10800886_2.fastq.gz,SRR10800886_3.fastq.gz

Save it into the :code:`Count_Basic/data` folder under :code:`experiment.csv`.

Running MPRAsnakeflow
---------------------

Now we have everything at hand to run the count MPRAsnakeflow pipeline. Therefore we have to be in the cloned MPRAsnakeflow folder. But we will change the working and output directory to the :code:`Count_Basic` folder in :code: `config.yaml` file. The MPRAsnakeflow count command is:


.. code-block:: bash

    cd <path/to/MPRAsnakeflow>/MPRAsnakeflow
    conda activate mprasnakeflow
    snakemake --configfile config/config.yaml --use-conda -p -c 4

.. note:: Please check your :code:`config/cluster.config` file and :code:`config/config.yaml` if it is correctly configured (e.g. with your SGE cluster commands).

If everything works fine the following 5 processes will run: :code:`create_BAM (make idx)` :code:`raw_counts`, :code:`filter_counts`, :code:`final_counts` and :code:`add-ons`, :code:`dna_rna_merge_counts`, :code:`calc_correlations`, :code:`make_master_tables`.


Results
-----------------

All output files will be in the :code:`results/(name of the project)` folder.

We expect the program to output the following status when complete:

.. code-block:: text

    Finished job 100.
    100 of 100 steps (100%) done

To generate a final report, the following code can be used

.. code-block:: bash

    snakemake --report report.html