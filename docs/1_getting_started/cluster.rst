.. _Cluster:

=======================================
Running MPRAsnakeflow on HPC Cluster
=======================================

Snakemake gives us the opportunity to run MPRAsnakeflow in a cluster environment. Please check the Snakemake documentation for more information on how to set up a cluster environment. We use `snakemake resources <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources>`_ to set the main resources per rule. Most resources are generic and can be used on multiple clusters, environments, or even locally. We have a predefined workflow profile with resources: :download:`config.yaml <../profiles/default/config.yaml>`:

.. literalinclude:: ../profiles/default/config.yaml
    :language: yaml
    :lines: 1-30

We used this workflow successfully in a SLURM environment using the `slurm executor plugin <https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html>`_ from Snakemake. Therefore, the partition is set with :code:`slurm_partition` and has to be renamed or removed to fit with your own SLURM configuration.

Running with resources
----------------------

Example: Using 30 cores and 10GB of memory.

.. code-block:: bash

    snakemake --sdm conda --configfile config/config.yaml -c 30 --resources mem_mb=10000 --workflow-profile profiles/default

Performance tweaks: Running specific rules with different resources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some rules will benefit from multithreading or more memory. This can be specified within your profile, workflow profile, or in the command line interface using :code:`--set-resources RULE_NAME:RESOURCE_NAME=VALUE` or :code:`--set-threads RULE_NAME=VALUE`. Before changing resources, make sure that you really need the rule by running a dry run to get the list of executed rules only:

.. code-block:: bash

    snakemake -n --quiet rules

Possible rules to tweak:

:Assignment:

    :assignment_hybridFWRead_get_reads_by_cutadapt:
        Only needed when using the linker option in the config. You can add more threads using :code:`--set-threads assignment_hybridFWRead_get_reads_by_cutadapt=4`. Default is always 1 thread.

    :assignment_mapping_bbmap:
        Only needed when using bbmap for mapping. Memory and threads can be optimized, e.g., via :code:`--set-threads assignment_mapping_bbmap=30 --set-resources assignment_mapping_bbmap:mem_mb=10000`. Default is 1 thread and 4GB memory, but we recommend using 30 threads and 10GB if available.

    :assignment_mapping_bwa:
        Only needed when using bwa for mapping. Memory and threads can be optimized, e.g., via :code:`--set-threads assignment_mapping_bwa=30 --set-resources assignment_mapping_bwa:mem_mb=10000`. Default is 1 thread, but we recommend using 30 threads and 10GB if available.

    :assignment_collectBCs:
        Threads can be optimized, e.g., via :code:`--set-threads assignment_collectBCs=30`. Default is 1 thread, but we recommend using 30 threads if available.

:Experiment:

    :counts_onlyFW_raw_counts_by_cutadapt:
        Only needed when you have only FW reads and use the adapter option. Threads can be optimized, e.g., via :code:`--set-threads experiment_counts_onlyFW_raw_counts_by_cutadapt=30`. Default is 1 thread.

Running on an HPC using SLURM
-----------------------------

Using the SLURM executor plugin to run 300 jobs in parallel:

.. code-block:: bash

    snakemake --sdm conda --configfile config/config.yaml -j 300 --workflow-profile profiles/default --executor slurm

Snakemake 7 (not supported anymore)
-----------------------------------

In Snakemake 7, we used the :code:`--cluster` option, which is not available in Snakemake 8. You can also use the predefined `config/sbatch.yaml`, but this might be outdated. We highly recommend using resources with the workflow profile.

.. code-block:: bash

    snakemake --use-conda --configfile config/config.yaml --cluster "sbatch --nodes=1 --ntasks={cluster.threads} --mem={cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output}" --jobs 100 --cluster-config config/sbatch.yaml

Please note that with this :code:`--cluster` option, the log folder of the cluster environment (see :code:`-o {cluster.output}`) has to be generated first, e.g.:

.. code-block:: bash

    mkdir -p logs

.. note:: Please consult your cluster's wiki page for cluster-specific commands and change cluster options to reflect these specifications. Additionally, for large libraries, more memory can be specified in this location.
