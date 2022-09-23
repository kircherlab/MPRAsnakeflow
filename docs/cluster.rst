.. _Cluster:

=======================================
Running MPRAsnakeflow on HPC Cluster
=======================================

Snakemake gives us the opportunity to run MPRAsnakeflow in a cluster environment. Please check the Snakemake documentation for more information on how to set up a cluster environment. We already configured a cluster environment for SLURM. This can be used to adapt the workflow to other HPC environments.

The configuration (SLURM) of resources for rules can be found in ``config/sbatch.yml``, allowing each process to be run as a separate ``sbatch`` command. For running MPRAsnakflow here is a possible command:

.. code-block:: bash

    snakemake --use-conda --configfile conf/config.yaml --cluster "sbatch --nodes=1 --ntasks={cluster.threads} --mem={cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output}" --jobs 100 --cluster-config config/sbatch.yaml

Please note that the log folder of the cluster environment has to be generated first, e.g:

.. code-block:: bash

    mkdir -p logs

  .. note:: Please consult your cluster's wiki page for cluster specific commands and change cluster Options to reflect these specifications. Additionally, for large libraries, more memory can be specified in this location.
