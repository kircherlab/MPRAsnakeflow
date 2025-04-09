.. _Installation:

=====================
Installation
=====================

Installation should take less than 5 minutes.

System Requirements
===================

- CentOS Linux 7 or above

Required Packages
=================

Package Management
------------------

You will need either `conda` or `apptainer` for managing dependencies:

1. **Conda** (version ≥ 24.7.1)

   .. code-block:: bash

       conda >24.7.1

   Download here: https://github.com/conda-forge/miniforge

2. **Apptainer** (for running software in containers)

   .. code-block:: bash

       apptainer

   Download here: https://apptainer.org/docs/user/latest/quick_start.html#installation

Workflow Language
-----------------

You will also need `snakemake` (version ≥ 8.24.1):

.. code-block:: bash

    snakemake 8.24.1

Download here: https://snakemake.readthedocs.io/

Clone Repository
=================

Clone the MPRAsnakeflow repository from GitHub:

.. code-block:: bash

    git clone https://github.com/kircherlab/MPRAsnakeflow.git

Set Up Snakemake Environment with Conda
=======================================

This pipeline uses Python 2.7, Python ≥3.7, and additional R scripts in a Snakemake pipeline. The provided `.yml` files will create the appropriate environments, which are completely handled by MPRAsnakeflow. The pipeline is designed to run on a Linux system.

To set up the environment:

1. Install the general conda environment called `snakemake`:

   .. code-block:: bash

       cd MPRAsnakeflow
       conda create -c conda-forge -c bioconda -n snakemake snakemake

2. Activate the `snakemake` environment:

   .. code-block:: bash

       conda activate snakemake

3. To deactivate the environment, use:

   .. code-block:: bash

       conda deactivate

Quick Test
==========

To verify the installation, run the following command:

.. code-block:: bash

    conda activate snakemake
    snakemake --help

