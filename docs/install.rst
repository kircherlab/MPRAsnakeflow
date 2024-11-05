.. _Installation:

=====================
Installation
=====================

Installation should take less than 5 minutes.

System Requirements
===================

CentOS Linux 7 or above

Required packages
=================

Package management
------------------

.. code-block:: bash

    conda >24.7.1 or above

Download here: https://github.com/conda-forge/miniforge

or for running software in containers

.. code-block:: bash

    apptainer

Download here: https://apptainer.org/docs/user/latest/quick_start.html#installation

Workflow language
-----------------

.. code-block:: bash

    snakemake 8.24.1 or above

Download here: https://snakemake.readthedocs.io/


Clone repository
=================

Download here: https://github.com/kircherlab/MPRAsnakeflow.git


Set up snakemake environment with conda
=============================================

This pipeline uses python2.7 and python ≥3.7 with additional R scripts in a Snakemake pipeline. The ``.yml`` files provided will create the appropriate environments and is completely handled by MPRAsnakeflow. The whole pipeline is set up to run on a Linux system.

Install the the conda environment. The general conda environment is called ``snakemake``.

.. code-block:: bash

    cd MPRAsnakeflow
    conda create -c conda-forge -c bioconda -n snakemake snakemake
    
    # activate snakemake
    conda activate snakemake

To deactivate the environment, use:

.. code-block:: bash

    conda deactivate



Quick test
============

.. code-block:: bash

    conda activate snakemake
    snakemake --help
    
