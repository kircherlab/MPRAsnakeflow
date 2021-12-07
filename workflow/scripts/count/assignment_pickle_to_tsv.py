# Author: Max Schubach 2021
# print out barcodes and correspoding assignment (tsv format to standard out) from a picke file

import pickle
import click

# options
@click.command()
@click.option('-i',
              '--input',
              'assignment_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Assignment pickle file.')
def cli(assignment_file):


    assoc=pickle.load(open(assignment_file,'rb'))


    for oligo,barcodes in assoc.items():
        for barcode in barcodes:
             print("%s\t%s" % (barcode, oligo))

if __name__ == '__main__':
    cli()
