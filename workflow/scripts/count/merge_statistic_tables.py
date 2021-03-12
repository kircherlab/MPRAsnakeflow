#Author Gracie Gordon 2019
#merge tables

import sys
import pandas as pd
import numpy as np

import click

# options
@click.command()
@click.option('--condition',
              required=True,
              type=str,
              help='Name of the condition.')
@click.option('--statistic',
              'statistic',
              required=True,
              nargs=2,
              multiple=True,
              type=(str, click.Path(exists=True, readable=True)),
              help='Replicate name and statistic file. Can be used multiple times')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file.')
def cli(condition, statistic, output_file):

    df_full = None
    for replicate_count in statistic:

        rep=replicate_count[0]
        file=replicate_count[1]

        df_stats=pd.read_csv(file, sep='\t', header=0)
        df_stats.insert(0,"condition",condition)
        df_stats.insert(1,"replicate",rep)

        if (df_full is not None):
            df_full=pd.concat([df_full,df_stats])
        else:
            df_full=df_stats

    df_full.to_csv(output_file, index=False,sep='\t', compression='gzip')


if __name__ == '__main__':
    cli()
