# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 10:16:29 2022

@author: dashpm

Samplerer- set fraction n so value is generated as follows.

usecase: samplerer.py --input <count_file> --prop <proportion to achieve> --output <output_file>
(accepted values as prop are float only between 0 and 1)
"""

import pandas as pd
import click
import math
import random

# options


@click.command()
@click.option('--input',
              ('input_file'),
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Barcodes count file (DNA or RNA) in .tsv format.')
@click.option('--prop',
              ('prop_val'),
              required=False, type=float,
              help='prop value e.g., 0.2, 0.3 (only between 0 and 1).')
@click.option('--total',
              ('total_val'),
              required=False, type=int,
              help='Total number of counts after sampling.')
@click.option('--threshold',
              ('threshold_val'),
              required=False, type=int,
              help='Set upper limit of counts allowed. e.g., 10, 20.')
@click.option('--seed',
              ('seed'),
              required=False, type=int,
              help='Use seed for random number generator')
@click.option('--minCounts',
              ('min_counts_val'),
              required=False, type=int,
              help='Remove all BCs with lower values.')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file.')
def cli(input_file, prop_val, total_val, threshold_val, seed, min_counts_val, output_file):
    # set seed if defined
    if seed is not None:
        random.seed(seed)
    # Filtering table
    click.echo("Reading count file...")
    df_ = pd.read_csv(input_file, header=None, sep='\t')
    if total_val is not None or prop_val is not None:
        # taking the smalles proportion when prop_val and total_val givebn
        pp = 1.0
        if prop_val is not None:
            pp = prop_val
        if total_val is not None:
            total_ = sum(df_.iloc[:, 1].values)
            pp = min(total_val/total_, pp)

        click.echo("Adjusting barcodes with given proportion %f" % pp)
        df_.iloc[:, 1] = df_.iloc[:, 1].astype(int).apply(lambda x: int(
            math.floor(x*pp) + (0.0 if random.random() > (x*pp-math.floor(x*pp)) else 1.0)))
    if threshold_val is not None:
        click.echo("Adjusting barcodes with counts > threshold...")
        df_.iloc[:, 1] = df_.iloc[:, 1].astype(int).apply(lambda x: threshold_val if x > threshold_val else x)

    if min_counts_val is not None:
        click.echo("Filter barcodes with min counts...")
        df_ = df_[df_.iloc[:, 1] >= min_counts_val]

    click.echo("Writing count file...")
    df_.to_csv(output_file, sep="\t", index=False, header=None, compression='gzip')


if __name__ == '__main__':
    cli()
