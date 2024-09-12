# Author: Max Schubach 2022

import pandas as pd
import numpy as np
import click


# options
@click.command()
@click.option('--counts',
              'counts_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Oligos with DNA and RNA counts as well as observed barcodes.')
@click.option('--declaration',
              'declaration_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help="Declaration file which is oligo is ref and alt. " +
              "Should be ID REF ALT as TSV (also with this header).")
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file with .')
def cli(counts_file, declaration_file, output_file):

    # declaration file
    click.echo("Read declaration file...")
    # expected column names: ID REF ALT
    declaration = pd.read_csv(declaration_file, header=0, sep="\t")

    # get count df
    click.echo("Read count file...")
    # expected column names: name dna_counts rna_counts dna_normalized rna_normalized ratio log2 n_obs_bc
    counts = pd.read_csv(counts_file, header=0, sep="\t", index_col=0)
    # join ref with index of counts
    output = declaration.join(counts, on='REF')
    # join again for the alt columns
    output = output.join(counts, on='ALT', lsuffix='_REF', rsuffix='_ALT')
    output["log2_expression"] = np.log2(output['ratio_ALT']/output['ratio_REF'])

    # fill NA and set correct output types
    output.fillna(0, inplace=True)
    output = output.astype(dtype={'dna_counts_REF': 'int64', 'rna_counts_REF': 'int64', 'n_obs_bc_REF': 'int64',
                                  'dna_counts_ALT': 'int64', 'rna_counts_ALT': 'int64', 'n_obs_bc_ALT': 'int64'},
                           copy=False)

    # write output
    click.echo("Write files...")
    output.to_csv(output_file, sep='\t', header=True, compression='gzip')


if __name__ == '__main__':
    cli()
