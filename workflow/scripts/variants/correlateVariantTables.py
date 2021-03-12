# Author: Max Schubach 20202

import sys
import pandas as pd
import numpy as np
import math
import click


# options
@click.command()
@click.option('--condition',
              required=True,
              type=str,
              help='Name of the condition.')
@click.option('--variants',
              'variants',
              required=True,
              nargs=2,
              multiple=True,
              type=(str, click.Path(exists=True, readable=True)),
              help='Replicate name and variat table file. Can be used multiple times')
@click.option('--bc-threshold',
              'bc_threshold',
              required=False,
              type=int,
              default=0,
              help='Number of observed barcodes thershold. Min. for ref and alt of both files')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file.')

def cli(condition, variants, bc_threshold, output_file):

    def filterOnThreshold(variants, threshold):
        return(variants.query('n_obs_bc_ALT >= %d & n_obs_bc_REF >= %d' % (threshold, threshold)))

    output = pd.DataFrame()
    for i in range(len(variants)):
        for j in range(i+1,len(variants)):

            rep_1=variants[i][0]
            file_1=variants[i][1]

            rep_2=variants[j][0]
            file_2=variants[j][1]

            click.echo("Compare replicate %s with replicate %s" % (rep_1, rep_2))

            # variants file
            click.echo("Read variants file...")
            variants_1 = pd.read_csv(file_1, header=0, sep="\t", index_col=0)
            variants_2 = pd.read_csv(file_2, header=0, sep="\t", index_col=0)

            if (bc_threshold > 0):
                click.echo("Filter variants file using min BC %d..." % bc_threshold)
                variants_1 = filterOnThreshold(variants_1, bc_threshold)
                variants_2 = filterOnThreshold(variants_2, bc_threshold)

            click.echo("Join variants file...")
            variants_join = variants_1.join(variants_2, how="inner", lsuffix='_A', rsuffix='_B')[["log2_expression_A", "log2_expression_B"]]




            output = output.append([[condition, rep_1, rep_2, variants_join.shape[0], bc_threshold, variants_join.corr(method="pearson").iloc[0,1],variants_join.corr(method="spearman").iloc[0,1]]])


    # write output
    click.echo("Write files...")
    output = output.rename(columns={0:"Condition", 1:"Replicate_A", 2:"Replicate_B", 3:"n_Variants", 4:"minBarcodes", 5:"Pearson", 6:"Spearman"})
    click.echo(output)
    output.to_csv(output_file, index=False, sep='\t', header=True, compression='gzip')


if __name__ == '__main__':
    cli()
