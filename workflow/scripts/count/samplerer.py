# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 10:16:29 2022

@author: dashpm

Samplerer- set mean and standard deviation so value is generated as follows.

usecase: samplerer.py --input <count_file> --mean <mean_value> --std <standard_deviation_value> --output <output_file>
(accepted values as mean/std are positive integers only)
"""

import pandas as pd
import click
import statistics as stat
import numpy as np

# options
@click.command()
@click.option('--input', 
              ('input_file'), 
              required=True, 
              type = click.Path(exists=True, readable=True), 
              help= 'Barcodes count file (DNA or RNA) in .tsv format.')
@click.option('--mean', 
              ('mean_val'), 
              required=True, type = int, 
              help= 'Mean value e.g., 10, 20.')
@click.option('--std', 
              ('std_val'), 
              required=True, type = int, 
              help= 'Standard deviation e.g., 5, 10.')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file.')
    
def cli(input_file, mean_val, std_val, output_file):
    
    # Filtering table 
    click.echo("Reading count file...")
    df_ = pd.read_csv(input_file, header=None, sep='\t')
    
    click.echo("Adjusting barcodes with given mean and standard deviation")
    dist_ = stat.NormalDist(mean_val, std_val)
    samples_ = dist_.samples(df_.shape[0], seed=42)  
    df_.iloc[:,1] = np.array(samples_, int)
    
    click.echo("Writing count file...")
    df_.to_csv(output_file, sep="\t", index=False, header=None, compression='gzip')
    
if __name__ == '__main__':
    cli()