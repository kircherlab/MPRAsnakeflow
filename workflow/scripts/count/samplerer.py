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
              type = click.Path(exists=True, readable=True), 
              help= 'Barcodes count file (DNA or RNA) in .tsv format.')
@click.option('--prop', 
              ('prop_val'), 
              required=False, type = float, 
              help= 'prop value e.g., 0.2, 0.3 (only between 0 and 1).')
@click.option('--threshold', 
              ('threshold_val'), 
              required=False, type = int, 
              help= 'Set upper limit of counts allowed. e.g., 10, 20.')
@click.option('--seed', 
              ('seed'), 
              required=False, type = int, 
              help= 'Use seed for random number generator')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file.')

def cli(input_file, prop_val, threshold_val, seed, output_file):
    # set seed if defined
    if seed:
        random.seed(seed)
    # Filtering table 
    click.echo("Reading count file...")
    df_ = pd.read_csv(input_file, header=None, sep='\t')
    if prop_val != None:
        total_ = sum(df_.iloc[:,1].values)
        pp = prop_val/total_
        click.echo("Adjusting barcodes with given proportion")
        df_.iloc[:,1] = df_.iloc[:,1].astype(int).apply(lambda x: int(math.floor(x*pp) + (0.0 if random.random() > (x*pp-math.floor(x*pp)) else 1.0)))
    if threshold_val != None:
        click.echo("Adjusting barcodes with counts > threshold...")
        df_.iloc[:,1] = df_.iloc[:,1].astype(int).apply(lambda x: threshold_val if x > threshold_val else x)   
    
    click.echo("Writing count file...")
    df_.to_csv(output_file, sep="\t", index=False, header=None, compression='gzip')
    
if __name__ == '__main__':
    cli()
