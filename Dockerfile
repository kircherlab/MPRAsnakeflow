ARG VERSION=0.5.0

FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="7a57714fe74eb25255d53b45e2095cd8a4dd4fe73db79006353670c432af97b1"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/NGmerge.yaml
#   prefix: /conda-envs/8b5bbd33d7ccdbbe3a28773771abe2b3
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - ngmerge=0.3
#     - python
#     - click
#     - htslib
RUN mkdir -p /conda-envs/8b5bbd33d7ccdbbe3a28773771abe2b3
COPY workflow/envs/NGmerge.yaml /conda-envs/8b5bbd33d7ccdbbe3a28773771abe2b3/environment.yaml

# Conda environment:
#   source: workflow/envs/bbmap_samtools_htslib.yaml
#   prefix: /conda-envs/f24b0ccfc23aadb93b466380cd592733
# ---
# channels:
#   - bioconda
#   - conda-forge
# dependencies:
#   - bbmap
#   - samtools
#   - htslib
RUN mkdir -p /conda-envs/f24b0ccfc23aadb93b466380cd592733
COPY workflow/envs/bbmap_samtools_htslib.yaml /conda-envs/f24b0ccfc23aadb93b466380cd592733/environment.yaml

# Conda environment:
#   source: workflow/envs/bwa_samtools_picard_htslib.yaml
#   prefix: /conda-envs/a0a27bbdca540d02023c7fdcf819bfc1
#   ---
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - bwa
#     - samtools
#     - picard
#     - htslib
RUN mkdir -p /conda-envs/a0a27bbdca540d02023c7fdcf819bfc1
COPY workflow/envs/bwa_samtools_picard_htslib.yaml /conda-envs/a0a27bbdca540d02023c7fdcf819bfc1/environment.yaml

# Conda environment:
#   source: workflow/envs/default.yaml
#   prefix: /conda-envs/899c0a8c04c6edeed4b242dbc3e0e1c8
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - htslib
RUN mkdir -p /conda-envs/899c0a8c04c6edeed4b242dbc3e0e1c8
COPY workflow/envs/default.yaml /conda-envs/899c0a8c04c6edeed4b242dbc3e0e1c8/environment.yaml

# Conda environment:
#   source: workflow/envs/fastqsplitter.yaml
#   prefix: /conda-envs/e5aec3a0d6b8921994e5305d4f89e90f
#   ---
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - fastqsplitter
RUN mkdir -p /conda-envs/e5aec3a0d6b8921994e5305d4f89e90f
COPY workflow/envs/fastqsplitter.yaml /conda-envs/e5aec3a0d6b8921994e5305d4f89e90f/environment.yaml

# Conda environment:
#   source: workflow/envs/python3.yaml
#   prefix: /conda-envs/a4e1b935cbca52df9b6f192ff86c464c
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - biopython
#     - click
#     - dask
#     - matplotlib
#     - numpy
#     - pandas
#     - python
#     - pysam
#     - pyfastx
RUN mkdir -p /conda-envs/a4e1b935cbca52df9b6f192ff86c464c
COPY workflow/envs/python3.yaml /conda-envs/a4e1b935cbca52df9b6f192ff86c464c/environment.yaml

# Conda environment:
#   source: workflow/envs/r.yaml
#   prefix: /conda-envs/ae3e37bf43cbb30416a885168e10c552
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - r-base
#     - r-cowplot
#     - r-cairo
#     - r-optparse
#     - r-tidyverse
RUN mkdir -p /conda-envs/ae3e37bf43cbb30416a885168e10c552
COPY workflow/envs/r.yaml /conda-envs/ae3e37bf43cbb30416a885168e10c552/environment.yaml

# Conda environment:
#   source: workflow/envs/python27.yaml
#   prefix: /conda-envs/cb972f023533b03e742da9095ce03b06
#   ---
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - htslib
#     - pysam
#     - python=2.7
#     - samtools
RUN mkdir -p /conda-envs/cb972f023533b03e742da9095ce03b06
COPY workflow/envs/python27.yaml /conda-envs/cb972f023533b03e742da9095ce03b06/environment.yaml

# Conda environment:
#   source: workflow/envs/cutadapt.yaml
#   prefix: /conda-envs/a3e2fce7f2f6fdbe1aa97232e3def601
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#    dependencies:
#     - cutadapt
RUN mkdir -p /conda-envs/a3e2fce7f2f6fdbe1aa97232e3def601
COPY workflow/envs/cutadapt.yaml /conda-envs/a3e2fce7f2f6fdbe1aa97232e3def601/environment.yaml

# Conda environment:
#   source: workflow/envs/quarto.yaml
#   prefix: /conda-envs/b933cc1aa7c25db04635e7ec0e37f80e
#   ---
#   channels:
#      - conda-forge
#      - bioconda
#   dependencies:
#      - python
#      - quarto
#      - jupyter
#      - pandas
#      - matplotlib
#      - papermill
RUN mkdir -p /conda-envs/b933cc1aa7c25db04635e7ec0e37f80e
COPY workflow/envs/quarto.yaml /conda-envs/b933cc1aa7c25db04635e7ec0e37f80e/environment.yaml

# Conda environment:
#   source: workflow/envs/mpralib.yaml
#   prefix: /conda-envs/f14db2daf3d94b49e39ea8fded7fe87e
#   ---
#   channels:
#       - bioconda
#       - conda-forge
#   dependencies:
#       - python
#       - click
#       - pip
#       - pandas
#       - numpy
#       - biopython
#       - matplotlib
#       - scikit-learn
#       - seaborn
#       - scipy
#       - anndata
#       - pysam
#       - scipy
#       - pip:
#           - mpralib==0.6.0
RUN mkdir -p /conda-envs/76616de101f1e6a8ce3a11b15175e85e
COPY workflow/envs/mpralib.yaml /conda-envs/76616de101f1e6a8ce3a11b15175e85e/environment.yaml


# Step 2: Generate conda environments

RUN <<EOR
	conda config --add channels nodefaults
	conda config --add channels bioconda
	conda config --add channels conda-forge
	conda config --set channel_priority strict
EOR

RUN conda env create --no-default-packages --prefix /conda-envs/8b5bbd33d7ccdbbe3a28773771abe2b3 --file /conda-envs/8b5bbd33d7ccdbbe3a28773771abe2b3/environment.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/f24b0ccfc23aadb93b466380cd592733 --file /conda-envs/f24b0ccfc23aadb93b466380cd592733/environment.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/a0a27bbdca540d02023c7fdcf819bfc1 --file /conda-envs/a0a27bbdca540d02023c7fdcf819bfc1/environment.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/a3e2fce7f2f6fdbe1aa97232e3def601 --file /conda-envs/a3e2fce7f2f6fdbe1aa97232e3def601/environment.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/899c0a8c04c6edeed4b242dbc3e0e1c8 --file /conda-envs/899c0a8c04c6edeed4b242dbc3e0e1c8/environment.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/e5aec3a0d6b8921994e5305d4f89e90f --file /conda-envs/e5aec3a0d6b8921994e5305d4f89e90f/environment.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/cb972f023533b03e742da9095ce03b06 --file /conda-envs/cb972f023533b03e742da9095ce03b06/environment.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/a4e1b935cbca52df9b6f192ff86c464c --file /conda-envs/a4e1b935cbca52df9b6f192ff86c464c/environment.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/b933cc1aa7c25db04635e7ec0e37f80e --file /conda-envs/b933cc1aa7c25db04635e7ec0e37f80e/environment.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/ae3e37bf43cbb30416a885168e10c552 --file /conda-envs/ae3e37bf43cbb30416a885168e10c552/environment.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/76616de101f1e6a8ce3a11b15175e85e --file /conda-envs/76616de101f1e6a8ce3a11b15175e85e/environment.yaml

# cleanup when version changed
ARG VERSION
RUN conda clean --all -y
