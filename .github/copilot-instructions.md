# MPRAsnakeflow Copilot Instructions

## Project Overview
MPRAsnakeflow is a Snakemake workflow (≥8.24.1) for processing Massively Parallel Reporter Assay (MPRA) sequencing data. It assigns barcodes to designed oligos and generates count tables.

## Architecture

### Core Components
- **Assignment workflow** ([workflow/rules/assignment.smk](workflow/rules/assignment.smk)): Maps barcodes to oligos using `exact`, `bwa`, or `bbmap` alignment tools
- **Experiment workflow** ([workflow/rules/experiment.smk](workflow/rules/experiment.smk)): Processes count data and merges with assignments
- **Common utilities** ([workflow/rules/common.smk](workflow/rules/common.smk)): Helper functions, config validation, wildcard constraints

### Data Flow
1. Config validation against JSON schema ([workflow/schemas/config.schema.yaml](workflow/schemas/config.schema.yaml))
2. Design file validation → FASTQ splitting → Alignment → Barcode collection → Assignment output → QC report
3. Experiment counts → Merge with assignments → QC report (Quarto-based)

### Key Wildcards
Constrained in Snakefile: `project`, `assignment`, `condition`, `replicate`, `type` (DNA/RNA), `config`

## Development Commands

```bash
# Lint Snakemake rules
snakemake --lint --configfile config/example_config.yaml --config skip_version_check=True

# Format Snakemake files
snakefmt .

# Dry run
snakemake --software-deployment-method conda --configfile config/example_config.yaml -n

# Build docs
cd docs && make html
```

## Code Conventions

### Snakemake Rules
- Use `getScript()` and `getCondaEnv()` helpers from [common.smk](workflow/rules/common.smk) for paths
- Each rule needs: `conda` directive (env from `workflow/envs/`), `input`, `output`, `log`, and `shell`/`script`
- Temp files: wrap with `temp()`, logs go to `results/logs/`

### Python Scripts
- Located in `workflow/scripts/` with subdirs per workflow component
- Use `click` for CLI arguments (see [check_design_file.py](workflow/scripts/assignment/check_design_file.py))
- Flake8 config: max-line-length=127, select E9/F6/F7/W4/W8

### Config Schema
- All config options must be defined in [config.schema.yaml](workflow/schemas/config.schema.yaml)
- Version field required; validated against MPRAsnakeflow version
- Conditional validation using JSON Schema `allOf`/`if`/`then` for tool-specific options

### Conda Environments
- One env file per tool group in `workflow/envs/`
- Use `mpralib.yaml` for mpralib-dependent scripts, `python3.yaml` for general Python

## CI/CD
- GitHub Actions: docs build, super-linter (Python/YAML/JSON/Snakefmt/R), Snakemake lint
- Testing is commented out but uses `snakemake-github-action`

## Common Patterns

### Adding a New Rule
```snakemake
rule my_rule:
    conda:
        getCondaEnv("python3.yaml")
    input:
        script=getScript("mydir/myscript.py"),
        data=lambda wc: config["assignments"][wc.assignment]["design_file"],
    output:
        "results/assignment/{assignment}/myoutput.tsv",
    log:
        "results/logs/assignment/my_rule.{assignment}.log",
    shell:
        "python {input.script} --input {input.data} > {output} 2> {log}"
```

### Config Access
Access via lambdas: `lambda wc: config["assignments"][wc.assignment]["alignment_tool"]["tool"]`
