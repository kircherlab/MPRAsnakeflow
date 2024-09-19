FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="7a57714fe74eb25255d53b45e2095cd8a4dd4fe73db79006353670c432af97b1"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/NGmerge.yaml
#   prefix: /conda-envs/c243bde7dc056785a077f6c33e56e8d6
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - ngmerge=0.3
#     - python
#     - click  
#     - htslib
RUN mkdir -p /conda-envs/c243bde7dc056785a077f6c33e56e8d6
COPY workflow/envs/NGmerge.yaml /conda-envs/c243bde7dc056785a077f6c33e56e8d6/environment.yaml

# Conda environment:
#   source: workflow/envs/bbmap_samtools_htslib.yaml
#   prefix: /conda-envs/575ebc82fb464fb2d0748323abbd3a13
# ---
# channels:
#   - bioconda
#   - conda-forge
#   - defaults
# dependencies:
#   - bbmap
#   - samtools
#   - htslib
RUN mkdir -p /conda-envs/575ebc82fb464fb2d0748323abbd3a13
COPY workflow/envs/bbmap_samtools_htslib.yaml /conda-envs/575ebc82fb464fb2d0748323abbd3a13/environment.yaml

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
#   source: workflow/envs/fastqsplitter.yaml
#   prefix: /conda-envs/dc242c7dafc90db387bc0290c31dc7ae
#   ---
#   channels:
#     - bioconda
#     - defaults
#     - conda-forge
#   dependencies:
#     - fastqsplitter
RUN mkdir -p /conda-envs/dc242c7dafc90db387bc0290c31dc7ae
COPY workflow/envs/fastqsplitter.yaml /conda-envs/dc242c7dafc90db387bc0290c31dc7ae/environment.yaml

# Conda environment:
#   source: workflow/envs/python3.yaml
#   prefix: /conda-envs/dadb883da8c83465d38f12e012df0cd0
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
#     - python
#     - pysam
#     - pyfastx
RUN mkdir -p /conda-envs/dadb883da8c83465d38f12e012df0cd0
COPY workflow/envs/python3.yaml /conda-envs/dadb883da8c83465d38f12e012df0cd0/environment.yaml

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

# Conda environment:
#   source: workflow/envs/python27.yaml
#   prefix: /conda-envs/c1d850971f4158052cd52615fbc1591a
#   ---
#   channels:
#     - bioconda
#     - defaults
#     - conda-forge
#   dependencies:
#     - htslib
#     - pysam
#     - python=2.7
#     - samtools
RUN mkdir -p /conda-envs/c1d850971f4158052cd52615fbc1591a
COPY workflow/envs/python27.yaml /conda-envs/c1d850971f4158052cd52615fbc1591a/environment.yaml

# Conda environment:
#   source: workflow/envs/cutadapt.yaml
#   prefix: /conda-envs/e6c048b22dbbbe081b8d18143c20afe3
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#    dependencies:
#     - cutadapt
RUN mkdir -p /conda-envs/d49adba2589cd2a66656b9298acdbece
COPY workflow/envs/cutadapt.yaml /conda-envs/d49adba2589cd2a66656b9298acdbece/environment.yaml

# Conda environment:
#   source: workflow/envs/quarto.yaml
#   prefix: /conda-envs/b8e51d222ab0d9caac2206a127729b1c
#   channels:
#      - conda-forge
#      - bioconda
#      - defaults
#   dependencies:
#      - python
#      - quarto
#      - jupyter
#      - pandas
#      - matplotlib
#      - papermill
RUN mkdir -p /conda-envs/b8e51d222ab0d9caac2206a127729b1c
COPY workflow/envs/quarto.yaml /conda-envs/b8e51d222ab0d9caac2206a127729b1c/environment.yaml



# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/c243bde7dc056785a077f6c33e56e8d6 --file /conda-envs/c243bde7dc056785a077f6c33e56e8d6/environment.yaml
RUN mamba env create --prefix /conda-envs/575ebc82fb464fb2d0748323abbd3a13 --file /conda-envs/575ebc82fb464fb2d0748323abbd3a13/environment.yaml
RUN mamba env create --prefix /conda-envs/f354d1f7a8fd64abb8ea8902ec91d399 --file /conda-envs/f354d1f7a8fd64abb8ea8902ec91d399/environment.yaml
RUN mamba env create --prefix /conda-envs/9444545a0ebc79ec516fa74514742720 --file /conda-envs/9444545a0ebc79ec516fa74514742720/environment.yaml
RUN mamba env create --prefix /conda-envs/dc242c7dafc90db387bc0290c31dc7ae --file /conda-envs/dc242c7dafc90db387bc0290c31dc7ae/environment.yaml
RUN mamba env create --prefix /conda-envs/dadb883da8c83465d38f12e012df0cd0 --file /conda-envs/dadb883da8c83465d38f12e012df0cd0/environment.yaml
RUN mamba env create --prefix /conda-envs/e6c048b22dbbbe081b8d18143c20afe3 --file /conda-envs/e6c048b22dbbbe081b8d18143c20afe3/environment.yaml
RUN mamba env create --prefix /conda-envs/c1d850971f4158052cd52615fbc1591a --file /conda-envs/c1d850971f4158052cd52615fbc1591a/environment.yaml
RUN mamba env create --prefix /conda-envs/d49adba2589cd2a66656b9298acdbece --file /conda-envs/d49adba2589cd2a66656b9298acdbece/environment.yaml
RUN mamba env create --prefix /conda-envs/b8e51d222ab0d9caac2206a127729b1c --file /conda-envs/b8e51d222ab0d9caac2206a127729b1c/environment.yaml
RUN mamba clean --all -y
