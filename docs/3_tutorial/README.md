# MPRAsnakeflow tutorial

This repository contains example data and python notebooks for running MPRAsnakeflow. You can run the notebooks locally or within colab.

## Introduction

MPRAsnakeflow is a pipeline that processes sequencing data from MPRA to create count tables for candidate sequences tested in the experiment.
The code can be found [here](https://github.com/kircherlab/MPRAsnakeflow) and an extensive documentation [here](https://mprasnakeflow.readthedocs.io).  
This readme explains the types of notebooks that are within the repository and type of example data as well as the structure of the input or MPRAsnakeflow.

If you have more questions or need help, please write an issue on the [MPRAsnakeflow_tutorial repository](https://github.com/kircherlab/MPRAsnakeflow_tutorial/issues). For general questions about MPRAsnakeflow, please write an issue on the [MPRAsnakeflow repository](https://github.com/kircherlab/MPRAsnakeflow/issues). You can also contact the main developer Max Schubach (<max.schubach@bih-charite.de>).


## Notebooks

- Assignment: `tutorial_assignment.ipynb`  
  This notebook is used to assign the barcodes to the sequences and to create a final file which maps barcodes to designed oligos, like `example_data/assignment/assignment_barcodes.default.tsv.gz`. The config file for the part of the workflow is `config_assignment.yaml`.
- Experiment: `tutorial_experiment.ipynb`  
  This notebook is used to generate activity measurements from DNA, barcode counts. The input file of the assignment workflow is needed to associate the barcodes to designed oligos. The config file for the part of the workflow is `config_experiment.yaml`.
- Combined: `tutorial_combined.ipynb`  
  This notebook runs both, the assignment and the experiment workflow. It does not use a pre-defined assignment file (maps barcodes to oligos). It uses directly the output of the assignment, configured in the config file. The advantage is that more jobs can be parallalized because the oligo/barcode map is needed later in the experiment workflow. Therfore the jobs before can already be run. This is usually helpfull if you run MPRAsnakeflow on a large HPC using a scheduling system like SLURM. The config file for the part of the workflow is `config_combined.yaml`.


## Example data

We designed a CNN-based deep neural network (DNN) to identify variants potentially affecting specific tissues or being tissue-agnostic across diverse human cell lines (e.g., HepG2, HEK293T, K562, WTC-11). We selected 120,000 variants based on high/low predictive effects from our DNN for MPRA experiments across all mentioned cell types.  

All candidate cis-regulatory sequences (cCREs) of length 200bp (flanked from both directions with a 15bp barcode) undergo MPRA in human hepatocellular carcinoma cell line (HepG2) and kidney epithelial (HEK293T) cells as 3 replicates. This example dataset contains ~1000 variants and ~100 negative control sequences tested in HepG2.

The example data is stored in the `example_data/` directory. The config files you need for MPRAsnakflow are within this folder called `config_assignment.yaml`, `config_experiment.yaml`, and `config_combined.yaml`. For more information on the config files see the [MPRAsnakeflow documentation](https://mprasnakeflow.readthedocs.io/en/latest/config.html).

## Running the tutorial notebooks 

### Local
- If you want to run the workflow locally with the example data it took on a unix system with a single core 20 minutes and required around 3GB of RAM.
- You need to install [Docker](https://www.docker.com) or [Apptainer](https://apptainer.org) to run the workflow locally.
- Also an installation of jupyter notebooks with libraries matplotlib and cv2 is needed.

### On google colab

Google Colab Notebooks are a free, cloud-based Jupyter notebook environment that allows you to write and execute Python code. A more indepth guide can be found [here](https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/what_is_jupyter.html). The notebooks are usually stored in your Google Drive and can be shared just like Google Docs or Sheets. You can mount your Google Drive to the notebook and access your files from there which will be covered in the following tutorial. But also you can run the notebooks directly from the github repository.

- Click on the "Open in Colab" button in the notebook you want to run.
- It should be ablte to run the notebook on a standard, free resource of colab.
