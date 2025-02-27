## General configurations
Create or adjust the `config/example_config.yaml` in the repository to your needs to configure the workflow execution. When running on a cluster environment you need a special [executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/), e.g. like [SLURM](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html), and use an adapted workflow profile (original `profiles/default/config.yaml`) to set the correct values (like slurm partitions).

Detailed information about the config file can be found in the [MPRAsnakeflow documentation](https://mprasnakeflow.readthedocs.io/latest/index.html) under [Config file](https://mprasnakeflow.readthedocs.io/latest/config.html).
