# Author: Max Schubach 2021
# print out barcodes and correspoding assignment (tsv format to standard out) from a picke file

import click
import csv
import pandas as pd


# options
@click.command()
@click.option('--experiment',
              'experiment_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Experiment file.')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(exists=True, writable=True),
              help='Output file.')
def cli(experiment_file, output_file):

    exp = pd.read_csv(experiment_file)

    exp_dna = pd.concat(
        [
            exp[["BC_DNA"]],
            exp[["Condition", "Replicate"]].agg("_".join, axis=1) + "_DNA",
        ],
        axis=1,
    ).rename(columns={"BC_DNA": "BC", 0: "Sample"})

    exp_rna = pd.concat(
        [
            exp[["BC_RNA"]],
            exp[["Condition", "Replicate"]].agg("_".join, axis=1) + "_RNA",
        ],
        axis=1,
    ).rename(columns={"BC_RNA": "BC", 0: "Sample"})

    exp_dna.append(exp_rna).to_csv(output_file, sep="\t", header=False, index=False)


if __name__ == '__main__':
    cli()
