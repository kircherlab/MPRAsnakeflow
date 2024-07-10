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

    bcs=set()
    oligos=set()
    other_bcs=set()
    ambiguous_bc=set()

    chunksize = 10 ** 6
    chunk = pd.read_csv(input_file, usecols=[0, 1], names=["BC", "Oligo"], sep="\t", chunksize=chunksize)
    for df in chunk:
        bcs.update(df.BC)
        oligos.update(df.Oligo)
        if 'other' in df.Oligo:
            other_bcs.update(df[df.Oligo == "other"].BC)
        if 'ambiguous' in df.Oligo:
            ambiguous_bc.update(df[df.Oligo == "ambiguous"].BC)

    n_bcs = len(bcs)

    n_other_bcs = len(other_bcs)
    
    n_ambiguous_bc = len(ambiguous_bc)

    n_assigned_bcs = n_bcs - n_other_bcs - n_ambiguous_bc


    n_matched_oligos = len(oligos)

    if "ambiguous" in oligos:
        n_matched_oligos -= 1
    if "other" in oligos:
        n_matched_oligos -= 1

    output = pd.DataFrame({"Counts": [n_bcs, n_assigned_bcs, n_matched_oligos, n_other_bcs, n_ambiguous_bc]}, index=[
                          "BCs", "Assigned BCs", "Assigned oligos", "Other BCs", "Ambiguous BCs"])

    output.to_csv(output_file, sep='\t', header=True, index=True)


if __name__ == '__main__':
    cli()
