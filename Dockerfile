FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="5421b5b4528c6825eb311803d5f114f63dcdc2de14fd6290d123a2c12d0c204b"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/bwa_samtools_picard_htslib.yaml
#   prefix: /conda-envs/f354d1f7a8fd64abb8ea8902ec91d399
#   ---
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - bwa
#     - samtools
#     - picard
#     - htslib
RUN mkdir -p /conda-envs/f354d1f7a8fd64abb8ea8902ec91d399
COPY workflow/envs/bwa_samtools_picard_htslib.yaml /conda-envs/f354d1f7a8fd64abb8ea8902ec91d399/environment.yaml

# Conda environment:
#   source: workflow/envs/default.yaml
#   prefix: /conda-envs/9444545a0ebc79ec516fa74514742720
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - htslib
RUN mkdir -p /conda-envs/9444545a0ebc79ec516fa74514742720
COPY workflow/envs/default.yaml /conda-envs/9444545a0ebc79ec516fa74514742720/environment.yaml

# Conda environment:
#   source: workflow/envs/python3.yaml
#   prefix: /conda-envs/26726def692c11bc66e0e8363b8e54c0
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - biopython
#     - click
#     - dask
#     - matplotlib
#     - numpy
#     - pandas
#     - polars
#     - python
#     - pysam
RUN mkdir -p /conda-envs/26726def692c11bc66e0e8363b8e54c0
COPY workflow/envs/python3.yaml /conda-envs/26726def692c11bc66e0e8363b8e54c0/environment.yaml

# Conda environment:
#   source: workflow/envs/r.yaml
#   prefix: /conda-envs/e6c048b22dbbbe081b8d18143c20afe3
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - r-base
#     - r-cowplot
#     - r-cairo
#     - r-optparse
#     - r-tidyverse
RUN mkdir -p /conda-envs/e6c048b22dbbbe081b8d18143c20afe3
COPY workflow/envs/r.yaml /conda-envs/e6c048b22dbbbe081b8d18143c20afe3/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/f354d1f7a8fd64abb8ea8902ec91d399 --file /conda-envs/f354d1f7a8fd64abb8ea8902ec91d399/environment.yaml && \
    mamba env create --prefix /conda-envs/9444545a0ebc79ec516fa74514742720 --file /conda-envs/9444545a0ebc79ec516fa74514742720/environment.yaml && \
    mamba env create --prefix /conda-envs/26726def692c11bc66e0e8363b8e54c0 --file /conda-envs/26726def692c11bc66e0e8363b8e54c0/environment.yaml && \
    mamba env create --prefix /conda-envs/e6c048b22dbbbe081b8d18143c20afe3 --file /conda-envs/e6c048b22dbbbe081b8d18143c20afe3/environment.yaml && \
    mamba clean --all -y
