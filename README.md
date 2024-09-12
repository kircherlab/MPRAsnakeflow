# Snakemake workflow: MPRAsnakeflow

[![Documentation Status](https://readthedocs.org/projects/mprasnakeflow/badge/?version=latest)](https://mprasnakeflow.readthedocs.io/en/latest/?badge=latest)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.2.1-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Tests](https://github.com/kircherlab/MPRAsnakeflow/actions/workflows/main.yml/badge.svg)](https://github.com/kircherlab/MPRAsnakeflow/actions/workflows/main.yml)

This pipeline processes sequencing data from Massively Parallel Reporter Assays (MPRA) to create count tables for candidate sequences tested in the experiment.

MPRAsnakeflow is built on top of [Snakemake](https://snakemake.readthedocs.io).

## Authors

* Max Schubach (@visze), Berlin Institute of Health at Charité -- Universitätsklinikum Berlin, [Computational Genome Biology Group](https://kircherlab.bihealth.org)

## Documentation

You can find an extensive documentation [here](https://mprasnakeflow.readthedocs.io)

## Tutorial

We have [tutorial](https://github.com/kircherlab/MPRAsnakeflow_tutorial) created in jupyter notebooks to run MPRAsnakeflow locally or within colab.

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above). Here is a very short description of the usage. Please look at the [documentation](https://mprasnakeflow.readthedocs.io) for more comprehensive usage. 

### Step 1: Obtain a copy of this workflow

[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local system, into the place where you have access on your compute nodes. It does not necessarily be the the same folder where you start your analysis, but it can.

### Step 2: Configure workflow
Create or adjust the `config/example_config.yaml` in the repository to your needs to configure the workflow execution. When running on a cluster environment you need a special [exccecutor plugin](https://snakemake.github.io/snakemake-plugin-catalog/), e.g. like [SLURM](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html), and use an adapted workflow profile (original `profiles/default/config.yaml`) to set the correct values (like slurm partitions).

### Step 3: Install Snakemake

Install Snakemake (recommended version >= 8.x) using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) (recommended installation via [miniforge](https://github.com/conda-forge/miniforge)):

    mamba create -c bioconda -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    mamba activate snakemake

Test your configuration by performing a dry-run via

    snakemake --software-deployment-method conda --configfile config.yaml -n

Execute the workflow locally via

    snakemake --software-deployment-method conda --cores $N --configfile config.yaml --workflow-profile profiles/default

using `$N` cores or run it in a cluster environment (here SLURM) via the [slurm excecutor plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html),

    snakemake --software-deployment-method conda --executor slurm --cores $N --configfile config.yaml --workflow-profile profiles/default

Please note that `profiles/default/config.yaml` has to be adapted to your needs (like partition names).
For snakemake 7.x this might work too using slurm sbatch (but depricated in newer snakemake versions:

    snakemake --use-conda --configfile config.yaml --cluster "sbatch --nodes=1 --ntasks={cluster.threads} --mem={cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output}" --jobs 100 --cluster-config config/sbatch.yaml


Please note that the log folder of the cluster environment has to be generated first, e.g:

    mkdir -p logs

For other cluster environments please check the [Snakemake](https://snakemake.readthedocs.io) documentation nad look for other [exccecutor plugins](https://snakemake.github.io/snakemake-plugin-catalog/) and adapt accodingly.

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --sdm apptainer,conda --cores $N --configfile config.yaml --workflow-profile profiles/default

in combination with any of the modes above. This will use a pre-build singularity container of MPRAsnakeflow with the conda ens installed in.


It is also possible to run the workflow in a different folder so that the results get stored not in the MPRAsnakeflow folder. Here you have to specify the snakefile path, like

    snakemake --sdm conda --configfile yourConfigFile.yaml --snakefile <path/to/MPRAsnakeflow>/MPRAsnakeflow/workflow/Snakefile --cores $N --workflow-profile <path/to/MPRAsnakeflow>/profiles/default

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

### Step 5: Investigate results

This part still works but it is outdated. Use the QC report, see documentation.
After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake --report report.html --configfile config.yaml

This report can, e.g., be forwarded to your collaborators.

### Step 6: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the repository:

    git commit -a
    git push

### Step 7: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: `git remote add -f upstream git@github.com:snakemake-workflows/MPRAsnakeflow.git` or `git remote add -f upstream https://github.com/snakemake-workflows/MPRAsnakeflow.git` if you do not have setup ssh keys.
2. Update the upstream version: `git fetch upstream`.
3. Create a diff with the current version: `git diff HEAD upstream/master workflow > upstream-changes.diff`.
4. Investigate the changes: `vim upstream-changes.diff`.
5. Apply the modified diff via: `git apply upstream-changes.diff`.
6. Carefully check whether you need to update the config files: `git diff HEAD upstream/master config`. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.


### Step 8: Contribute back

In case you have also changed or added steps, please consider contributing them back to the original repository:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the original repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to your local system, to a different place than where you ran your analysis.
3. Copy the modified files from your analysis to the clone of your fork, e.g., `cp -r workflow path/to/fork`. Make sure to **not** accidentally copy config file contents or sample sheets. Instead, manually update the example config files if necessary.
4. Commit and push your changes to your fork.
5. Create a [pull request](https://help.github.com/en/articles/creating-a-pull-request) against the original repository.

## Testing

Test cases are in the subfolder `.test`. They are automatically executed via continuous integration with [Github Actions](https://github.com/features/actions).

