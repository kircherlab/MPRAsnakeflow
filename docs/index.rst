.. _Homepage:

====================================
MPRAsnakeflow's documentation
====================================

.. image:: https://img.shields.io/badge/snakemake-≥8.24.1-brightgreen.svg
    :target: https://snakemake.github.io/

.. image:: https://img.shields.io/badge/conda->24.7.1-brightgreen.svg
    :target: https://github.com/conda-forge/miniforge


**Welcome!**

MPRAsnakeflow pipeline processes sequencing data from Massively Parallel Reporter Assays (MPRAs)
to create count tables for candidate sequences tested in the experiment.

MPRAsnakeflow is built on top of `Snakemake <https://snakemake.readthedocs.io/>`_ (version ≥8.24.1 required) and is configured via a ``.yaml`` file.

Authors
    Max Schubach (`@visze <https://github.com/visze>`_)
    `Computational Genome Biology Group <https://kircherlab.bihealth.org>`_
    Berlin Institute of Health at Charité - Universitätsklinikum Berlin 

Usage
    If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of the (original) repository.

Installation & Getting Started
    Instructions for the Installation of the program and some examples to get you started.

MPRAsnakeflow Workflows
    An overview of how MPRAsnakeflow works and documentation for the MPRAsnakeflow sub-workflows.

MPRAsnakeflow Turorials
    Get to know MPRAsnakeflow by runnig it via jupyter notebooks or colab on small examples.

MPRAsnakeflow Examples
    Multiple examples from the literature are listed for every sub-workflow in MPRAsnakeflow.
    
Tips & Tricks
    Find our FAQ here.

Project Information
    More information on the project, including the changelog, list of contributing authors, and contribution instructions.


-------------
Quick Example
-------------

To run MPRAsnakeflow, first activate the snakemake 8 environment with the following command:

.. code-block:: bash

    conda activate snakemake


And then run the main workflow with:

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
     - When ``conda`` is set, the utility uses conda to efficiently query repositories and query package dependencies. MPRAsnakeflow also can use containers via apptainer by using ``--software-deployment-method apptainer conda``. This will use a container to run all rules but inside it will activate the pre-installed conda environments. Recommended option: ``--software-deployment-method apptainer conda``
   * - ``--cores``
     - This utility sets the number of cores (``$N``) to be used by MPRAsnakeflow.
   * - ``--configfile``
     - This file (e.g., ``config/example_config.yaml``) contains the project, its objects and properties, and sub-properties and its objects that **must** be set before running MPRAsnakeflow.


-------------------
Investigate results
-------------------

Th ebest option to investiate your results is to have a look at the QC report.

(In development) After successful execution, you can create a self-contained interactive HTML report with all results via:

.. code-block:: bash

    snakemake --report report.html --configfile conf/example_config.yaml


This report can be forwarded to your collaborators.
An example of a generated report (using some trivial test data) can be seen `here <https://snakemake.github.io/resources/report.html>`_.

--------
Feedback
--------

Feel free to leave feedback(s), ask question(s), or report bug(s) at our issues page: `MPRAsnakeflow Issues <https://github.com/kircherlab/MPRAsnakeflow/issues>`_.


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

.. toctree::
    :caption: Installation & Getting Started
    :name: getting-started
    :maxdepth: 1
    :hidden:

    quickstart
    install
    config
    cluster


.. toctree::
    :caption: MPRAsnakeflow Workflows
    :name: mprasnakeflow-workflows
    :maxdepth: 1
    :hidden:

    overview
    assignment
    experiment

.. toctree::
    :caption: MPRAsnakeflow Tutorials
    :name: mprasnakeflow-tutorials
    :maxdepth: 1
    :hidden:

    tutorial

.. toctree::
    :caption: MPRAsnakeflow Examples
    :name: mprasnakeflow-examples
    :maxdepth: 2
    :hidden:

    assignment_example1
    count_example1
    combined_example1


.. toctree::
    :caption: Tips & Tricks
    :name: tips-tricks
    :maxdepth: 1
    :hidden:

    faq



.. toctree::
   :caption: Project Info
   :name: project-info
   :maxdepth: 1
   :hidden:

   contributing
   authors
   history
   license
   todo_list
