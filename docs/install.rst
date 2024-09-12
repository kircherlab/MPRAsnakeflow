.. _Installation:

=====================
Installation
=====================

Installation should take less than 5 minutes

System Requirements
===================

CentOS Linux 7 or above

Required packages
=================

Package management
------------------

.. code-block:: bash

    conda (mamba) 4.6 or above

Download here: https://github.com/conda-forge/miniforge

or for running software in containers

.. code-block:: bash

    apptainer

Download here: https://apptainer.org/docs/user/latest/quick_start.html#installation

Workflow language
-----------------

.. code-block:: bash

    snakemake 8.16.0 or above (snakemake >=7.15.1 will also work but cli might be different as here documented)

Download here: https://snakemake.readthedocs.io/


Clone repository
=================

Download here: https://github.com/kircherlab/MPRAsnakeflow.git


Set up snakemake environment with conda/mamba
=============================================

This pipeline uses python2.7 and python3.6 with additional R scripts in a Snakemake pipeline. The ``.yml`` files provided will create the appropriate environments and is completely handled by MPRAsnakeflow. The whole pipeline is set up to run on a Linux system.

Install the the conda environment. The general conda environment is called ``snakemake``.

.. code-block:: bash

    cd MPRAsnakeflow
    mamba create -c conda-forge -c bioconda -n snakemake snakemake
    
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
    
