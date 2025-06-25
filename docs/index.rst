.. _Homepage:

====================================
MPRAsnakeflow's Documentation
====================================

.. image:: https://img.shields.io/badge/snakemake-≥8.24.1-brightgreen.svg
    :target: https://snakemake.github.io/

.. image:: https://img.shields.io/badge/conda->24.7.1-brightgreen.svg
    :target: https://github.com/conda-forge/miniforge

**Welcome!**

MPRAsnakeflow is a pipeline designed to process sequencing data from Massively Parallel Reporter Assays (MPRAs) to create count tables for candidate sequences tested in the experiment.

MPRAsnakeflow is built on top of `Snakemake <https://snakemake.readthedocs.io/>`_ (version ≥8.24.1 required) and is configured via a ``.yaml`` file.

Authors
-------
- Max Schubach (`@visze <https://github.com/visze>`_) | `Computational Genome Biology Group <https://kircherlab.bihealth.org>`_ | Berlin Institute of Health at Charité - Universitätsklinikum Berlin 

Usage
-----
If you use this workflow in a paper, don't forget to give credit to the authors by citing the URL of the (original) repository.

Installation & Getting Started
------------------------------
Instructions for installing the program and examples to help you get started.

MPRAsnakeflow Workflows
-----------------------
An overview of how MPRAsnakeflow works and documentation for the MPRAsnakeflow sub-workflows.

MPRAsnakeflow Tutorials
-----------------------
Get to know MPRAsnakeflow by running it via Jupyter notebooks or Google Colab on small examples.

MPRAsnakeflow Examples
----------------------
Multiple examples from the literature are listed for every sub-workflow in MPRAsnakeflow.

Tips & Tricks
-------------
Find our FAQ here.

Project Information
-------------------
More information on the project, including the changelog, list of contributing authors, and contribution instructions.

-------------
Quick Example
-------------

To run MPRAsnakeflow, first activate the Snakemake environment with the following command:

.. code-block:: bash

    conda activate snakemake

Then run the main workflow with:

.. code-block:: bash

    snakemake --software-deployment-method conda --cores $N --configfile config/example_config.yaml

--------
Features
--------

.. list-table:: 
   :widths: 25 80
   :header-rows: 1

   * - Option
     - Description
   * - ``--software-deployment-method``
     - When ``conda`` is set, the utility uses conda to efficiently query repositories and manage package dependencies. MPRAsnakeflow can also use containers via Apptainer by using ``--software-deployment-method apptainer conda``. This will use a container to run all rules but will activate the pre-installed conda environments inside it. Recommended option: ``--software-deployment-method apptainer conda``.
   * - ``--cores``
     - This utility sets the number of cores (``$N``) to be used by MPRAsnakeflow.
   * - ``--configfile``
     - This file (e.g., ``config/example_config.yaml``) contains the project configuration, including objects and properties that **must** be set before running MPRAsnakeflow.

-------------------
Investigate Results
-------------------

The best option to investigate your results is to review the QC report.

(In development) After successful execution, you can create a self-contained interactive HTML report with all results via:

.. code-block:: bash

    snakemake --report report.html --configfile config/example_config.yaml

This report can be shared with collaborators. An example of a generated report (using some trivial test data) can be seen `here <https://snakemake.github.io/resources/report.html>`_.

--------
Feedback
--------

Feel free to leave feedback, ask questions, or report bugs on our issues page: `MPRAsnakeflow Issues <https://github.com/kircherlab/MPRAsnakeflow/issues>`_.

Indices and Tables
==================

* :ref:`genindex`
* :ref:`search`

.. toctree::
    :caption: Installation & Getting Started
    :name: getting-started
    :maxdepth: 1
    :hidden:

    1_getting_started/quickstart
    1_getting_started/install
    1_getting_started/config
    1_getting_started/cluster

.. toctree::
    :caption: MPRAsnakeflow Workflows
    :name: mprasnakeflow-workflows
    :maxdepth: 1
    :hidden:

    2_workflows/overview
    2_workflows/assignment
    2_workflows/experiment

.. toctree::
    :caption: MPRAsnakeflow Tutorials
    :name: mprasnakeflow-tutorials
    :maxdepth: 1
    :hidden:

    3_tutorial/tutorial

.. toctree::
    :caption: MPRAsnakeflow Examples
    :name: mprasnakeflow-examples
    :maxdepth: 2
    :hidden:

    4_examples/assignment_example1
    4_examples/count_example1
    4_examples/combined_example1
    4_examples/plasmid_example

.. toctree::
    :caption: Tips & Tricks
    :name: tips-tricks
    :maxdepth: 1
    :hidden:

    5_tips/faq

.. toctree::
   :caption: Project Info
   :name: project-info
   :maxdepth: 1
   :hidden:

   6_info/contributing
   6_info/authors
   6_info/history
   6_info/license
   6_info/todo_list
