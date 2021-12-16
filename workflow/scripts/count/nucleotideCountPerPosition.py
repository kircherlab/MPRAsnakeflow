# Author: Max Schubach 2021
# Nucleotides count per position

import pandas as pd
import numpy as np

import click


# options
@click.command()
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="TSV file with a column of nucleotides.",
)
@click.option(
    "--column", "column", required=True, type=str, help="Column name of number"
)
@click.option(
    "--header/--no-header",
    "header",
    default=False,
    show_default=True,
    help="If set the file contains header names.",
)
@click.option(
    "--output",
    "output_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file.",
)
def cli(input_file, column, header, output_file):
    if header:
        df = pd.read_csv(input_file, sep="\t")
        df = pd.DataFrame(df[column])
        df.rename(columns={column : "DNA"}, inplace=True)
    else:
        df = pd.read_csv(input_file, header=None, sep="\t")
        df = pd.DataFrame(df[int(column)-1])
        df.rename(columns={(int(column)-1) : "DNA"}, inplace=True)

    nucleotides = df["DNA"].apply(lambda x: pd.Series(list(x)))

    output = pd.DataFrame(index=["A","T","C","G","N"])

    count_columns = []
    for i, columns in nucleotides.iteritems():
        output["Position_%d_counts" % (i+1)] = nucleotides[i].value_counts()
        count_columns += ["Position_%d_counts" % (i+1)]
        output["Position_%d_ratio" % (i+1)] = nucleotides[i].value_counts(normalize=True)

    output.fillna(0, inplace=True)
    for column in count_columns:
        output[column] = pd.to_numeric(output[column], downcast='integer')
    output.index.name = "Nucleotide"
    output.to_csv(output_file, index=True, sep='\t')


if __name__ == "__main__":
    cli()
