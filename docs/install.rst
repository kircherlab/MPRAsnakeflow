.. _Installation:

=====================
Installation
=====================

Installation should take less than 5 minutes

System Requirements
===================

CentOS Linux 7 or above

Required packages
==================

.. code-block:: bash

  	conda (mamba) 4.6 or above

Download here: https://docs.conda.io/en/latest/miniconda.html

.. code-block:: bash

    snakemake 7.15.1 or above

Download here: https://snakemake.readthedocs.io/


Clone repository
=================

Download here: https://github.com/kircherlab/MPRAsnakeflow.git


Set up conda with snakemake environment
==========================================

This pipeline uses python2.7 and python3.6 with additional R scripts in a Snakemake pipeline. The ``.yml`` files provided will create the appropriate environments and is completely handled by MPRAsnakeflow. The whole pipeline is set up to run on a Linux system.

Install the the conda environment. The general conda environment is called ``mprasnakeflow``.

.. code-block:: bash

    cd MPRAsnakeflow
    conda create -c bioconda -c conda-forge -n mprasnakeflow snakemake
    
    # activate MPRAsnakeflow
    conda activate mprasnakeflow

To deactivate the environment, use:

.. code-block:: bash

    conda deactivate



Quick test
============

.. code-block:: bash

    conda activate mprasnakeflow
    snakemake --help
    
