.. _Cluster:

=======================================
Running MPRAsnakeflow on HPC Cluster
=======================================

Snakemake gives us the opportunity to run MPRAsnakeflow in a cluster environment. Please check the Snakemake documentation for more information on how to set up a cluster environment. We use `snakemake resources <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources>`_ to set the main resources per rule. Most resources are generic and can be used on multipe clusters, environments or even local. We have a preddefined workflow profile with resources: :download:`config.yaml <../profiles/default/config.yaml>`:

.. literalinclude:: ../workflow/schemas/config.schema.yaml
    :language: yaml
    :lines: 1-30


We used this workflow successfully in a SLURM environment using the `slurm excecutor plugin <https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html>`_ from snakemake. therfore the partition is set with :code:`slurm_partition` and has to be renamed maybe due to your environment.

Running with resources
----------------------

having 30 cores and 10GB of memory.

.. code-block:: bash

    snakemake --sdm conda --configfile config/config.yaml -c 30 --resources mem_mb=10000  --workflow-profile profiles/default --executor slurm

Running on an HPC using SLURM
-----------------------------

Using the slurm excecutor plugin running 300 jobs in parallel.

.. code-block:: bash

    snakemake --sdm conda --configfile config/config.yaml -j 300  --workflow-profile profiles/default --executor slurm


Snakemake 7
===========

Here we used the :code:`cluster` option which is not anymore avialable in snakemake 8. You can also use the predefined `config/sbatch.yaml` but this might be outdated and we highly recommend to use resources with the workfloe profile. 

.. code-block:: bash

    snakemake --use-conda --configfile config/config.yaml --cluster "sbatch --nodes=1 --ntasks={cluster.threads} --mem={cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output}" --jobs 100 --cluster-config config/sbatch.yaml

Please note that the log folder of the cluster environment has to be generated first, e.g:

.. code-block:: bash

    mkdir -p logs

.. note:: Please consult your cluster's wiki page for cluster specific commands and change cluster Options to reflect these specifications. Additionally, for large libraries, more memory can be specified in this location.
