#Author: Gracie Gordon 2019
#merge tables

import pandas as pd
import click


# options
@click.command()
@click.option('--counts',
              'counts_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Barcodes (DNA or RNA).')
@click.option('--assignment',
              'assignment_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Assignment tsv file.')
@click.option('--name',
              'name',
              required=True,
              type=str,
              help='Name the statistic should be written with.')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file.')
@click.option('--statistic',
              'statistic_file',
              required=True,
              type=click.Path(writable=True),
              help='Statistic output file.')
def cli(counts_file, assignment_file, output_file, statistic_file, name):


    # statistic
    statistic = pd.DataFrame(data={'Experiment' : [name], 'Barcodes': [0], 'Counts' : [0], 'Average counts' : [0],
                    'Assigned barcodes' : [0], 'Assigned counts' : [0], 'Average assigned counts' : [0],
                    'Fraction assigned barcodes' : [0], 'Fraction assigned counts' : [0]})
    # Association file
    click.echo("Read assignment file...")
    assoc_barcodes_oligo=pd.read_csv(assignment_file, header=None, usecols=[0,1], sep="\t", names=['Barcode','Oligo'])
    assoc_barcodes = set(assoc_barcodes_oligo.Barcode)


    #get count df
    click.echo("Read count file...")
    counts=pd.read_csv(counts_file, header=None, sep="\t", names=['Barcode','Counts'])

    statistic['Barcodes'] = counts.shape[0]
    statistic['Counts'] = sum(counts.Counts)
    statistic['Average counts'] = statistic['Counts']/statistic['Barcodes']

    #fill in labels from dictionary
    click.echo("Collecting rows to remove...")
    remove_idx=[]
    for i,row in counts.iterrows():
        if row.Barcode not in assoc_barcodes:
            remove_idx.append(i)

    click.echo("Remove not assigned...")
    counts.drop(remove_idx,inplace=True)

    statistic['Assigned barcodes'] = counts.shape[0]
    statistic['Assigned counts'] = sum(counts.Counts)
    statistic['Average assigned counts'] = statistic['Assigned counts']/statistic['Assigned barcodes']

    statistic['Fraction assigned barcodes'] = statistic['Assigned barcodes']/statistic['Barcodes']
    statistic['Fraction assigned counts'] = statistic['Assigned counts']/statistic['Counts']

    click.echo("Write files...")
    counts.to_csv(output_file, index=False, sep='\t', header=False, compression='gzip')

    statistic.to_csv(statistic_file, index=False,sep='\t', compression='gzip')

if __name__ == '__main__':
    cli()
