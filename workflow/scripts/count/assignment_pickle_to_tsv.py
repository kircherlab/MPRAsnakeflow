import pandas as pd
import pickle


#Author: Gracie Gordon 2019
#merge tables

import sys
import pandas as pd
import numpy as np

import math
import pickle

from Bio import SeqIO

import click


# options
@click.command()
@click.option('--input',
              'assignment_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Assignment pickle file.')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file.')
def cli(assignment_file, output_file):


    assoc=pickle.load(open(assignment_file,'rb'))

    statistic['oligos design'] = len(assoc)

    BC_key = {}
    for k,v in assoc.items():
        statistic[['barcodes design']] += len(v)
        for x in v:
             BC_key.setdefault(x,k)

if __name__ == '__main__':
    cli()
