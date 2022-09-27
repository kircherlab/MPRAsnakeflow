.. _Assignment example:

.. role:: bash(code)
   :language: bash

============================
Basic assignment workflow
============================

This example runs the assignment workflow on 5'/5' WT MRPA data in the HEPG2 cell line from `Klein J., Agarwal, V., Keith, A., et al. 2019 <https://www.biorxiv.org/content/10.1101/576405v1.full.pdf>`_.

Prerequirements
======================

This example depends on the following data and software:


Installation of MPRAsnakeflow
----------------------------------------

Please install conda, the MPRAsnakeflow environment and clone the actual MPRAsnakeflow master branch. You will find more help under :ref:`Installation`.

Meta Data
___________

It is necessary to get the ordered oligo array so that each enhancer sequence can be labeled in the analysis and to trim any adaptors still in the sequence, in this case we trim off 15bp from the end of each sequence

.. code-block:: bash

    mkdir -p Assoc_Basic/data
    cd Assoc_Basic/data
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4237nnn/GSM4237954/suppl/GSM4237954_9MPRA_elements.fa.gz

    zcat GSM4237954_9MPRA_elements.fa.gz |awk '{ count+=1; if (count == 1) { print } else { print substr($1,1,171)}; if (count == 2) { count=0 } }' > design.fa

Reads
----------

There is one set of association sequencing for this data, which contains a forward (CRS-forward), reverse (CRS-reverse), and index (barcode) read for DNA and RNA. These data must be downloaded. All data is publically available on the short read archive (SRA). We will use SRA-toolkit to obtain the data.

.. note:: You need 10 GB disk space to download the data!

.. code-block:: bash

    conda install sra-tools
    cd Assoc_Basic/data
    fastq-dump --gzip --split-files SRR10800986
    cd ..

For large files and unstable internet connection we reccommend the comand `prefetch` from SRA tools before running `fastq-dump`. This command is much smarter in warnings when something went wrong.

.. code-block:: bash

    conda install sra-tools
    cd Assoc_Basic/data
    prefetch SRR10800986
    fastq-dump --gzip --split-files SRR10800986
    cd ..

.. note:: Please be sure that all files are downloaded completely without errors! Depending on your internet connection this can take a while. If you just want some data to run MPRsnakeAflow you can just limit yourself to one condition and/or just one replicate. 


With

.. code-block:: bash

    tree data


the folder should look like this:

.. code-block:: text

    data

Here is an overview of the files:

.. csv-table:: HEPG2 association data
   :header: "Condition", "GEO Accession", "SRA Accession", SRA Runs
   :widths: 40, 10, 10, 20

   "HEPG2-association: HEPG2 library association", GSM4237954, SRX7474872, "SRR10800986"


MPRAsnakeflow
=================================

Now we are ready to run MPRAsnakeflow and create CRS-barcode mappings.

Run snakemake
------------------------------

Now we have everything at hand to run the count MPRAsnakeflow pipeline. We will run the pipeline directly in the :code:`Assoc_Basic` folder. The MPRAsnakeflow workflow can be in a different directory. Let's assume :code:`/home/user/MPRAsnakeflow`. 

First we have to configure the config file:

.. todo:: config file of assignment example

The MPRAsnakeflow command is:


.. code-block:: bash

    cd Assoc_Basic
    conda activate mprasnakeflow
    snakemake -c 1 --use-conda --snakefile /home/user/MPRAsnakeflow/workflow/Snakefile --config config.yml

.. note:: Please modify your code when running in a cluster environment. We have an example SLURM config file here :code:`config/sbatch.yml`.

If everything works fine the following 7 rules will run: :code:`count_bc_nolab` :code:`create_BWA_ref`, :code:`PE_merge`, :code:`align_BWA_PE`, :code:`collect_chunks`, :code:`map_element_barcodes`, :code:`filter_barcodes`.

.. todo:: Rules not correct in example assignment workflow


Results
-----------------

All needed output files will be in the :code:`results/assignemnts/assoc_basic` folder.
