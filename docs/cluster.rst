.. _Cluster:

=======================================
Running MPRAsnakeflow on HPC Cluster
=======================================

Snakemake gives us the opportunity to run MPRAsnakeflow in a cluster environment. Right now we split up processes into two main groups: ``longtime`` and ``shorttime``. We can define different job setting for both groups. As you can imagine from the names `longtime` defines processes that takes a while when running. Sometimes several days. `shortime` defines processes that are quicker and are usually done in several minutes.

To enable the submission to your cluster you have to edit the ``config/sbatch.yml``, allowing each process to be run as a separate ``sbatch`` or command. If you run MPRAsnakeflow on a cluster system make sure be that you use the ``sbatch_cluster_snakemake`` command and adapt to the options by following the ``--help`` help command.

  .. note:: Please consult your cluster's wiki page for cluster specific commands and change cluster Options = to reflect these specifications. Additionally, for large libraries, more memory can be specified in this location.
