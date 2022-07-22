# Author: Max Schubach 2022

import pandas as pd
import click


# options
@click.command()
@click.option('--input',
              'input_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='file with BC oligo/other counts')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file.')
def cli(input_file, output_file):

    df = pd.read_csv(input_file, usecols=[0, 1], names=["BC", "Oligo"], sep="\t")

    bcs = df.BC.unique().size

    oligos = df.Oligo.unique()

    if 'other' in oligos:
        other_bcs = df[df.Oligo == "other"].BC.unique().size
    else:
        other_bcs = 0

    if 'ambiguous' in oligos:
        ambiguous_bc = df[df.Oligo == "ambiguous"].BC.unique().size
    else:
        ambiguous_bc = 0

    df = df[df.Oligo != "other"]
    df = df[df.Oligo != "ambiguous"]

    assigned_bcs = df.BC.unique().size

    matched_oligos = df.Oligo.unique().size

    output = pd.DataFrame({"Counts": [bcs, assigned_bcs, matched_oligos, other_bcs, ambiguous_bc]}, index=[
                          "Total BCs", "Total assigned BCs", "Total assigned oligos", "Total other BCs", "Total ambiguous BCs"])

    output.to_csv(output_file, sep='\t', header=True, index=True)


if __name__ == '__main__':
    cli()
