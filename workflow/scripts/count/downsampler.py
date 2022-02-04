# -*- coding: utf-8 -*-
"""
Created on Thu Jan 9 10:16:29 2022

@author: dashpm

Thresholder - set a threshold and all counts above max_threshold is filtered

usecase: thresholder.py <count_file> <threshold_value>
(accepted values as threshold are positive integers only)
"""

import pandas as pd
import click

# options
@click.command()
@click.option('--input', 
              ('input_file'), 
              required=True, 
              type = click.Path(exists=True, readable=True), 
              help= 'Barcodes count file (DNA or RNA) in .tsv format.')
@click.option('--threshold', 
              ('threshold_val'), 
              required=True, type = int, 
              help= 'Maximum number of barcodes allowed. e.g., 10, 20.')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file.')

def cli(input_file, threshold_val, output_file):
    
    # Filtering table 
    click.echo("Reading count file...")
    df_ = pd.read_csv(input_file, header=None, sep='\t')
    
    click.echo("Adjusting barcodes with counts > threshold...")
    df_.iloc[:,1] = df_.iloc[:,1].astype(int).apply(lambda x: threshold_val if x > threshold_val else x)
    
    click.echo("Writing count file...")
    df_.to_csv(output_file, sep="\t", index=False, header=None, compression='gzip')
    
if __name__ == '__main__':
    cli()