ARG VERSION=0.6.0

FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="9a2cd7d943a12a67aae93978ccc10c6213b846480366ddb973a181835e0b7947"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/NGmerge.yaml
#   prefix: /conda-envs/2abb0048c6dce1e9bf1a7960f16f3a74
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - ngmerge=0.3
#     - python
#     - click
#     - htslib==1.21
RUN mkdir -p /conda-envs/2abb0048c6dce1e9bf1a7960f16f3a74
COPY workflow/envs/NGmerge.yaml /conda-envs/2abb0048c6dce1e9bf1a7960f16f3a74/environment.yaml

# Conda environment:
#   source: workflow/envs/fastq-join.yaml
#   prefix: /conda-envs/7f3db13e2aa951f4484dd79393a5b358
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - fastq-join=1.3.1
RUN mkdir -p /conda-envs/7f3db13e2aa951f4484dd79393a5b358
COPY workflow/envs/fastq-join.yaml /conda-envs/7f3db13e2aa951f4484dd79393a5b358/environment.yaml

# Conda environment:
#   source: workflow/envs/bbmap_samtools_htslib.yaml
#   prefix: /conda-envs/a1cc34886525015a2366c351dc84f094
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - bbmap
#     - samtools==1.21
#     - htslib==1.21
RUN mkdir -p /conda-envs/a1cc34886525015a2366c351dc84f094
COPY workflow/envs/bbmap_samtools_htslib.yaml /conda-envs/a1cc34886525015a2366c351dc84f094/environment.yaml

# Conda environment:
#   source: workflow/envs/bwa_samtools_picard_htslib.yaml
#   prefix: /conda-envs/bd140014e61f38f39c52be56455899ed
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - bwa
#     - samtools==1.21
#     - picard
#     - htslib==1.21
RUN mkdir -p /conda-envs/bd140014e61f38f39c52be56455899ed
COPY workflow/envs/bwa_samtools_picard_htslib.yaml /conda-envs/bd140014e61f38f39c52be56455899ed/environment.yaml

# Conda environment:
#   source: workflow/envs/default.yaml
#   prefix: /conda-envs/bc0b5d56a68abe252a6b7f205988f848
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - htslib==1.21
RUN mkdir -p /conda-envs/bc0b5d56a68abe252a6b7f205988f848
COPY workflow/envs/default.yaml /conda-envs/bc0b5d56a68abe252a6b7f205988f848/environment.yaml

# Conda environment:
#   source: workflow/envs/fastqsplitter.yaml
#   prefix: /conda-envs/e9dc603fe82e22dfa110fb43a9265aca
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - fastqsplitter
RUN mkdir -p /conda-envs/e9dc603fe82e22dfa110fb43a9265aca
COPY workflow/envs/fastqsplitter.yaml /conda-envs/e9dc603fe82e22dfa110fb43a9265aca/environment.yaml

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
#   prefix: /conda-envs/05a8eaa56b4a44e5531fbad1610b05d8
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - htslib
#     - pysam
#     - python=2.7
#     - samtools
RUN mkdir -p /conda-envs/05a8eaa56b4a44e5531fbad1610b05d8
COPY workflow/envs/python27.yaml /conda-envs/05a8eaa56b4a44e5531fbad1610b05d8/environment.yaml

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
#   prefix: /conda-envs/8a065417982df56be3b310bca28e7765
#   ---
#   channels:
#       - bioconda
#       - conda-forge
#   dependencies:
#       - 'conda-forge::biopython>=1.80'
#       - mpralib=0.10.3
RUN mkdir -p /conda-envs/8a065417982df56be3b310bca28e7765
COPY workflow/envs/mpralib.yaml /conda-envs/8a065417982df56be3b310bca28e7765/environment.yaml

# Conda environment:
#   source: workflow/envs/pbmm2_pysam.yaml
#   prefix: /conda-envs/2308b21c334f9613fdb840777a17d2b9
#   ---
#   channels:
#       - conda-forge
#       - bioconda
#   dependencies:
#       - pbmm2
#       - pysam
#       - biopython
#       - python>=3.10
RUN mkdir -p /conda-envs/2308b21c334f9613fdb840777a17d2b9
COPY workflow/envs/pbmm2_pysam.yaml /conda-envs/2308b21c334f9613fdb840777a17d2b9/environment.yaml

# Step 2: Generate conda environments

RUN <<EOR
	conda config --add channels nodefaults
	conda config --add channels bioconda
	conda config --add channels conda-forge
	conda config --set channel_priority strict
EOR

# workflow/envs/NGmerge.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/2abb0048c6dce1e9bf1a7960f16f3a74 --file /conda-envs/2abb0048c6dce1e9bf1a7960f16f3a74/environment.yaml
# workflow/envs/fastq-join.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/7f3db13e2aa951f4484dd79393a5b358 --file /conda-envs/7f3db13e2aa951f4484dd79393a5b358/environment.yaml
# workflow/envs/bbmap_samtools_htslib.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/a1cc34886525015a2366c351dc84f094 --file /conda-envs/a1cc34886525015a2366c351dc84f094/environment.yaml
# workflow/envs/bwa_samtools_picard_htslib.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/bd140014e61f38f39c52be56455899ed --file /conda-envs/bd140014e61f38f39c52be56455899ed/environment.yaml
# workflow/envs/cutadapt.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/a3e2fce7f2f6fdbe1aa97232e3def601 --file /conda-envs/a3e2fce7f2f6fdbe1aa97232e3def601/environment.yaml
# workflow/envs/default.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/bc0b5d56a68abe252a6b7f205988f848 --file /conda-envs/bc0b5d56a68abe252a6b7f205988f848/environment.yaml
# workflow/envs/fastqsplitter.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/e9dc603fe82e22dfa110fb43a9265aca --file /conda-envs/e9dc603fe82e22dfa110fb43a9265aca/environment.yaml
# workflow/envs/python27.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/05a8eaa56b4a44e5531fbad1610b05d8 --file /conda-envs/05a8eaa56b4a44e5531fbad1610b05d8/environment.yaml
# workflow/envs/python3.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/a4e1b935cbca52df9b6f192ff86c464c --file /conda-envs/a4e1b935cbca52df9b6f192ff86c464c/environment.yaml
# workflow/envs/quarto.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/b933cc1aa7c25db04635e7ec0e37f80e --file /conda-envs/b933cc1aa7c25db04635e7ec0e37f80e/environment.yaml
# workflow/envs/r.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/ae3e37bf43cbb30416a885168e10c552 --file /conda-envs/ae3e37bf43cbb30416a885168e10c552/environment.yaml
# workflow/envs/mpralib.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/8a065417982df56be3b310bca28e7765 --file /conda-envs/8a065417982df56be3b310bca28e7765/environment.yaml
# workflow/envs/pbmm2_pysam.yaml
RUN conda env create --no-default-packages --prefix /conda-envs/2308b21c334f9613fdb840777a17d2b9 --file /conda-envs/2308b21c334f9613fdb840777a17d2b9/environment.yaml

# cleanup when version changed
ARG VERSION
RUN conda clean --all -y
