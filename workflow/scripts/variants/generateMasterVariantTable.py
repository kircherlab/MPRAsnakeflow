# Author: Max Schubach 2022

import pandas as pd
import numpy as np
import click


# options
@click.command()
@click.option('--input',
              'input_files',
              multiple=True,
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Variant file for each replicate.')
@click.option('--minRNACounts',
              'minRNACounts',
              default=1,
              show_default=True,
              type=int,
              help='Minimum number of RNA counts required per barcode. If 0 pesudocounts are used.')
@click.option('--minDNACounts',
              'minDNACounts',
              default=0,
              show_default=True,
              type=int,
              help='Minimum number of DNA counts required per barcode. If 0 pesudocounts are used.')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file with combined variant expression (of all replicates).')
def cli(input_files, minRNACounts, minDNACounts, output_file):

    pseudocountDNA = 1 if minDNACounts == 0 else 0
    pseudocountRNA = 1 if minRNACounts == 0 else 0

    df = pd.DataFrame()
    # variant files
    click.echo("Read variant files per replicate...")
    for input_file in input_files:
        df = pd.concat([df, pd.read_csv(input_file, header=0, sep="\t")], axis=0, sort=False)

    # filter based on minDNACounts and minRNACounts
    df = df.query('dna_counts_REF >= %d & dna_counts_ALT >= %d' % (minDNACounts, minDNACounts))

    df = df.query('rna_counts_REF >= %d & rna_counts_ALT >= %d' % (minRNACounts, minRNACounts))

    click.echo("Create new expression values...")
    df = df.groupby(['ID', 'REF', 'ALT']).agg(dna_counts_REF=('dna_counts_REF', sum),
                                              rna_counts_REF=('rna_counts_REF', sum),
                                              n_bc_REF=('n_bc_REF', sum),
                                              dna_counts_ALT=('dna_counts_ALT', sum),
                                              rna_counts_ALT=('rna_counts_ALT', sum),
                                              n_bc_ALT=('n_bc_ALT', sum)
                                              ).reset_index()

    scaling = 10**6

    df_total = df[['dna_counts_REF', 'rna_counts_REF', 'n_bc_REF',
                   'dna_counts_ALT', 'rna_counts_ALT', 'n_bc_ALT']].sum()

    total_bc = df_total['n_bc_REF'] + df_total['n_bc_ALT']
    total_dna = df_total['dna_counts_REF'] + df_total['dna_counts_ALT'] + total_bc * pseudocountDNA
    total_rna = df_total['rna_counts_REF'] + df_total['rna_counts_ALT'] + total_bc * pseudocountRNA

    scaling = 10**6

    def normalize(df, total, count_type, ref_alt, pseudocount, scaling):
        return(((
            df["%s_counts_%s" % (count_type, ref_alt)] + pseudocount * df['n_bc_%s' % ref_alt]
        )/df['n_bc_%s' % ref_alt])/total*scaling)

    for ref_alt in ["REF", "ALT"]:
        for count_type in ["dna", "rna"]:
            total = total_rna if count_type == 'rna' else total_dna
            pseudocount = pseudocountRNA if count_type == 'rna' else pseudocountDNA

            df['%s_normalized_%s' % (count_type, ref_alt)] = normalize(
                df, total, count_type, ref_alt, pseudocount, scaling)

        df['ratio_%s' % ref_alt] = df['rna_normalized_%s' % ref_alt] / df['dna_normalized_%s' % ref_alt]
        df['log2FoldChange_%s' % ref_alt] = np.log2(df['ratio_%s' % ref_alt])

    df["log2FoldChange_expression"] = np.log2(df['ratio_ALT']/df['ratio_REF'])

    # fill NA and set correct output types
    df.fillna(0, inplace=True)
    df = df.astype(dtype={'dna_counts_REF': 'int64', 'rna_counts_REF': 'int64', 'n_bc_REF': 'int64',
                          'dna_counts_ALT': 'int64', 'rna_counts_ALT': 'int64', 'n_bc_ALT': 'int64'}, copy=False)

    df = df.reindex(columns=["ID", "REF", "ALT", "dna_counts_REF", "rna_counts_REF",
                             "dna_normalized_REF", "rna_normalized_REF", "ratio_REF", "log2FoldChange_REF",
                             "n_bc_REF", "dna_counts_ALT", "rna_counts_ALT",
                             "dna_normalized_ALT", "rna_normalized_ALT", "ratio_ALT",
                             "log2FoldChange_ALT", "n_bc_ALT", "log2FoldChange_expression"])

    # write output
    click.echo("Write files...")
    df.to_csv(output_file, index=False, sep='\t', header=True, compression='gzip')


if __name__ == '__main__':
    cli()
