#contributed by Tal Ashuach, Max Schubach

import os
import re
import sys
import gzip
import pandas as pd
import numpy as np
import click
from collections import defaultdict

# options
@click.command()
@click.option('--input',
              'input_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Merged DNA and RNA counts of all replicates.')
@click.option('--rna-counts-output',
              'rna_counts_output_file',
              required=True,
              type=click.Path(writable=True),
              help='RNA counts output for MPRAnalyze')
@click.option('--dna-counts-output',
              'dna_counts_output_file',
              required=True,
              type=click.Path(writable=True),
              help='DNA counts output for MPRAnalyze.')
@click.option('--rna-annotation-output',
              'rna_annotation_output_file',
              required=True,
              type=click.Path(writable=True),
              help='RNA annotation output for MPRAnalyze.')
@click.option('--dna-annotation-output',
              'dna_annotation_output_file',
              required=True,
              type=click.Path(writable=True),
              help='DNA annotation output for MPRAnalyze.')



def cli(input_file, rna_counts_output_file, dna_counts_output_file, rna_annotation_output_file, dna_annotation_output_file):

    annot_pattern = re.compile("^([DR]NA).*\(condition (.*), replicate (.*)\)$")
    def get_annot(head):
        m = annot_pattern.match(head)
        if m is not None:
            return m.group(1,2,3)

    # read input
    df = pd.read_csv(input_file,sep="\t", header='infer')

    # remove unknown BC
    df = df[df['label']!='no_BC'].set_index('label').fillna(0)


    # get replicates and DNA/RNAs
    header_annot = [get_annot(h) for h in list(df.columns)[2:]]
    dna_annot = pd.DataFrame([h for h in header_annot if h[0] == 'DNA']).rename(columns={0: "type", 1: "condition", 2 : "replicate"})
    rna_annot = pd.DataFrame([h for h in header_annot if h[0] == 'RNA']).rename(columns={0: "type", 1: "condition", 2 : "replicate"})
    n_dna_obs = len(dna_annot)
    n_rna_obs = len(rna_annot)

    # counts for observation

    dna_df = df.iloc[:,2:(2+n_dna_obs)].applymap(np.int64)
    rna_df = df.iloc[:,(2+n_dna_obs):].applymap(np.int64)

    ## generate output DNA/RNA annotations (type_condition_replicate_barcode)
    n_bc = df.groupby('label').Barcode.agg(len).max()
    def generateAnnotationOutput(data, number_barcodes):
        data = data.loc[data.index.repeat(number_barcodes)]
        data['barcode'] = data.groupby(['type','condition','replicate']).cumcount() +1
        data['barcode'] = data['barcode'].astype(str)
        data['sample'] = data[['type','condition','replicate','barcode']].agg('_'.join,axis=1)
        data = data[["sample", "type", "condition", "replicate", "barcode"]]
        return(data)
    dna_annot = generateAnnotationOutput(dna_annot, n_bc)
    rna_annot = generateAnnotationOutput(rna_annot, n_bc)

    ## generate output DNA/RNA count tables
    ## rows oligo/seq ids,/assignment then per barcode the counts. padding with zeros
    def generateCountOutput(data,columns):
        counts = pd.DataFrame(list(data.groupby('label').apply(lambda x: x.values.flatten()))).fillna(0).astype(np.int64)
        counts.columns = columns
        counts['seq_id'] = data.index.unique()
        counts = counts[(['seq_id'] + list(columns))]
        return(counts)

    dna_counts = generateCountOutput(dna_df,dna_annot['sample'])
    rna_counts = generateCountOutput(rna_df,rna_annot['sample'])

    ## write table function
    def write(data,file):
        data.to_csv(file, index=False,sep='\t', compression='gzip')

    ## write output DNA/RNA annotations
    write(dna_annot,dna_annotation_output_file)
    write(rna_annot,rna_annotation_output_file)

    ## write output DNA/RNA annotations
    write(dna_counts,dna_counts_output_file)
    write(rna_counts,rna_counts_output_file)

if __name__ == '__main__':
    cli()
