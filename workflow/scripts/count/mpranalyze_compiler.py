#contributed by Tal Ashuach, Max Schubach

import re

import click
import numpy as np
import pandas as pd

ANNOT_PATTERN = re.compile(r"^([DR]NA).*\(condition (.*), replicate (.*)\)$")


def get_annot(head: str) -> tuple[str|None, str|None, str|None]:
    match = ANNOT_PATTERN.match(head)
    if match is not None:
        group1 = match.group(1)
        group2 = match.group(2)
        group3 = match.group(3)
        return (group1, group2, group3)
    return (None, None, None)


def generate_annotation_output(data, number_barcodes):
    data = data.loc[data.index.repeat(number_barcodes)].copy()
    data['barcode'] = data.groupby(['type', 'condition', 'replicate']).cumcount() + 1
    data['barcode'] = data['barcode'].astype(str)
    data['sample'] = data[['type', 'condition', 'replicate', 'barcode']].agg('_'.join, axis=1)
    return data[["sample", "type", "condition", "replicate", "barcode"]]


def generate_count_output(data, columns, number_barcodes):
    rows = []
    seq_ids = []
    for label, group in data.groupby('label', sort=False):
        padded = np.zeros((number_barcodes, data.shape[1]), dtype=np.int64)
        vals = group.values[:number_barcodes].astype(np.int64)
        padded[:len(vals)] = vals
        rows.append(padded.flatten(order='F'))
        seq_ids.append(label)
    counts = pd.DataFrame(rows, columns=columns)
    counts.insert(0, 'seq_id', seq_ids)
    return counts


def write_table(data, file):
    data.to_csv(file, index=False, sep='\t', compression='gzip')


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

    dna_df = df.iloc[:,2:(2+n_dna_obs)].astype(np.int64)
    rna_df = df.iloc[:,(2+n_dna_obs):].astype(np.int64)

    ## generate output DNA/RNA annotations (type_condition_replicate_barcode)
    n_bc = df.groupby('label').Barcode.agg(len).max()
    dna_annot = generate_annotation_output(dna_annot, n_bc)
    rna_annot = generate_annotation_output(rna_annot, n_bc)

    ## generate output DNA/RNA count tables
    ## rows oligo/seq ids,/assignment then per barcode the counts. padding with zeros
    dna_counts = generate_count_output(dna_df, dna_annot['sample'], n_bc)
    rna_counts = generate_count_output(rna_df, rna_annot['sample'], n_bc)

    ## write output DNA/RNA annotations
    write_table(dna_annot, dna_annotation_output_file)
    write_table(rna_annot, rna_annotation_output_file)

    ## write output DNA/RNA annotations
    write_table(dna_counts, dna_counts_output_file)
    write_table(rna_counts, rna_counts_output_file)

if __name__ == '__main__':
    cli()
