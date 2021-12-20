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
    "use_header",
    default=False,
    show_default=True,
    help="If set the file contains header names.",
)
@click.option(
    "--chunksize",
    "chunksize",
    required=False,
    default=10000,
    type=int,
    help="Chunksize whiel reading.",
)
@click.option(
    "--output",
    "output_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file.",
)
def cli(input_file, column, use_header, output_file, chunksize):

    output = pd.DataFrame(index=["A", "T", "C", "G", "N"])
    count_columns = []
    ratio_columns = []

    header = 'infer' if use_header else None
    with pd.read_csv(input_file, header=header, sep="\t", chunksize=chunksize) as reader:
        for chunk in reader:

            if use_header:
                df = pd.DataFrame(chunk[column])
                df.rename(columns={column: "DNA"}, inplace=True)
            else:
                df = pd.DataFrame(chunk[int(column)-1])
                df.rename(columns={(int(column)-1): "DNA"}, inplace=True)

            nucleotides = df["DNA"].apply(lambda x: pd.Series(list(x)))

            for i, columns in nucleotides.iteritems():
                counts_column = "Position_%d_counts" % (i+1)
                ratio_column = "Position_%d_ratio" % (i+1)
                if counts_column not in output.columns:
                    output[counts_column] = 0.0
                    output[ratio_column] = 0.0
                    count_columns += [counts_column]
                    ratio_columns += [ratio_column]

                output[counts_column] = output[counts_column] + \
                    pd.Series(nucleotides[i].value_counts(), index=["A", "T", "C", "G", "N"]).fillna(0.0)

    for i in range(len(ratio_columns)):
        output[ratio_columns[i]] = output[count_columns[i]]/output[count_columns[i]].sum()
    for column in count_columns:
        output[column] = pd.to_numeric(output[column], downcast='integer')
    output.index.name = "Nucleotide"
    output.to_csv(output_file, index=True, sep='\t')


if __name__ == "__main__":
    cli()
