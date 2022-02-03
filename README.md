# Snakemake workflow: MPRAsnakeflow

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.32.0-brightgreen.svg)](https://snakemake.bitbucket.io)
![Tests](https://github.com/visze/MPRAsnakeflow/workflows/Tests/badge.svg)

[![Build Status](https://travis-ci.org/snakemake-workflows/MPRAsnakeflow.svg?branch=master)](https://travis-ci.org/snakemake-workflows/MPRAsnakeflow)

This pipeline processes sequencing data from Massively Parallel Reporter Assays (MPRA) to create count tables for candidate sequences tested in the experiment.

This is the template for a new Snakemake workflow. Replace this text with a comprehensive description covering the purpose and domain.
Insert your code into the respective folders, i.e. `scripts`, `rules`, and `envs`. Define the entry point of the workflow in the `Snakefile` and the main configuration in a `config.yaml` file.

## Authors

* Max Schubach (@visze), Berlin Institute of Health at Charité -- Universitätsklinikum Berlin, [Computational Genome Biology Group](https://kircherlab.bihealth.org)

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).

### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Create or adjust the `config.yaml` to configure the workflow execution. When running on a cluster environment there are `drmaa.yaml` for drmaa runs or `cluster.yaml` for SLURM environment which contain the resources required for each job.

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda --configfile conf/config.yaml -n

Execute the workflow locally via

    snakemake --use-conda --cores $N --configfile conf/config.yaml

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --configfile conf/config.yaml --cluster "sbatch --nodes=1 --ntasks={cluster.threads} --mem={cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output}" --jobs 100 --cluster-config config/cluster.yaml

or

    snakemake --use-conda --configfile conf/config.yaml --drmaa "-n {cluster.threads} --mem={cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output}" --jobs 100
    
Please note that the log folder of the cluster environment has to be generated first, e.g:

    mkdir -p logs

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity --configfile conf/config.yaml

in combination with any of the modes above.

It is also possible to run the workflow in a different folder so that the results get stored not in the MPRAsnakeflow folder. Here you have to specify the snakefile path, like

    snakemake --use-conda --configfile yourConfigFile.yaml --snakefile <path/to/MPRAsnakeflow>/MPRAsnakeflow/workflow/Snakefile --cores $N

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

### Step 5: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake --report report.html

This report can, e.g., be forwarded to your collaborators.
An example (using some trivial test data) can be seen [here](https://cdn.rawgit.com/snakemake-workflows/rna-seq-kallisto-sleuth/master/.test/report.html).

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

