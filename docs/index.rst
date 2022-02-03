====================================
Welcome to MPRAsnakeflow's documentation!
====================================

.. image:: https://img.shields.io/badge/snakemake-≥5.32.0-brightgreen.svg
    :target: https://snakemake.bitbucket.io


This pipeline processes sequencing data from Massively Parallel Reporter Assays (MPRA)
to create count tables for candidate sequences tested in the experiment.

This is the template for a new Snakemake workflow. Replace this text with a comprehensive description covering the purpose and domain. Insert your code into the respective folders, i.e. ``scripts``, ``rules``, and ``envs``. Define the entry point of the workflow in the ``Snakefile`` and the main configuration in a ``config.yaml`` file.

Authors
    Max Schubach (`@visze <https://github.com/visze>`_), Berlin Institute of Health at Charité -- Universitätsklinikum Berlin, `Computational Genome Biology Group <https://kircherlab.bihealth.org>`_.

Usage
    If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of the (original) repository and, if available, its DOI (see above).

Installation & Getting Started
    Instructions for the Installation of the program and some examples to get you started.

MPRAsnakeflow Workflows
    An overview of how MPRAsnakeflow works and documentation for the MPRAsnakeflow sub workflows.

MPRAsnakeflow Examples
    Muliple examples from the literature are listed for every sub workflow in MPRAsnakeflow.

Project Information
    More information on the project, including the changelog, list of contributing authors, and contribution instructions.


-------------
Quick Example
-------------

.. code-block:: bash

    snakemake --use-conda --cores $N --configfile config/config.yaml


--------
Features
--------

:--use-conda:
  This utility uses mamba to efficiently query repositories and query package dependencies.
:--cores:
  This utility sets the number of cores to be used by MPRAsnakeflow.
:--configfile:
  This file (e.g., ``config/config.yaml``) contains the project, its objects and properties, and sub-properties and its objects that must be set before running MPRAsnakeflow.


Investigate results
    After successful execution, you can create a self-contained interactive HTML report with all results via:

.. code-block:: bash

    snakemake --report report.html


This report can, e.g., be forwarded to your collaborators.
An example (using some trivial test data) can be seen `here <https://cdn.rawgit.com/snakemake-workflows/rna-seq-kallisto-sleuth/master/.test/report.html>`_.

--------
Feedback
--------

The best place to leave feedback, ask questions, and report bugs is the `MPRAsnakeflow Issue Tracker <https://github.com/kircherlab/MPRAsnakeflow/issues>`_.


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
    cluster


.. toctree::
    :caption: MPRAsnakeflow Workflows
    :name: mprasnakeflow-workflows
    :maxdepth: 1
    :hidden:

    overview
    association
    count
    association_saturation_mutagenesis
    saturation_mutagenesis

.. toctree::
    :caption: MPRAsnakeflow Examples
    :name: mprasnakeflow-examples
    :maxdepth: 2
    :hidden:

    association_example1
    count_example1
    count_example2
    association_saturation_mutagenesis_example1
    saturation_mutagenesis_example1


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
