# -*- coding: utf-8 -*-
"""
Created on Wed Apr 06 14:00:09 2022

@author: Max Schubach

Samplerer - set fraction n so value is generated as follows.

usecase: samplerer.py --input <assignment_file> --prop <proportion to achieve> --output <output_file>
"""

import pandas as pd
import click

# options


@click.command()
@click.option('--input',
              ('input_file'),
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Assignment file in .tsv format.')
@click.option('--prop',
              ('prop_val'),
              required=False, type=float,
              help='Prop value e.g., 0.2, 0.3 (only between 0 and 1).')
@click.option('--total',
              ('total_val'),
              required=False, type=int,
              help='Total number of barcodes after sampling.')
@click.option('--seed',
              ('seed'),
              required=False, type=int,
              help='Use seed for random number generator')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file.')
def cli(input_file, prop_val, total_val, seed, output_file):
    # Filtering table
    click.echo("Reading assignment file...")
    df_ = pd.read_csv(input_file, header=None, sep='\t')

    click.echo("Removing duplicate assignments...")
    df_.drop_duplicates(subset=[0], keep = False, inplace=True)

    if total_val is not None or prop_val is not None:
        # taking the smalles proportion when prop_val and total_val givebn
        pp = 1.0
        if prop_val is not None:
            pp = prop_val
        if total_val is not None:
            total_ = len(df_)
            pp = min(total_val/total_, pp)
        
        click.echo("Adjusting barcodes with given proportion %f" % pp)
        df_ = df_.sample(frac = pp, replace = False, random_state = seed).sort_values(by=[1,0])

    click.echo("Writing count file...")
    df_.to_csv(output_file, sep="\t", index=False, header=None, compression='gzip')


if __name__ == '__main__':
    cli()
