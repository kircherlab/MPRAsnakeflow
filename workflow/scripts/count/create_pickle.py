###
# python create_pickle.py /fast/groups/ag_kircher/MPRA/MPRA_ML_Experiments/assignment/assignment_barcodes_incl_other.sorted.tsv.gz assignment_barcodes_incl_other.sorted.pickle
###


import click
import pickle
import gzip
import csv

@click.command()
@click.option("-i",
              "--input",
              'input_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help="Assignment tsv file (barcode TAB oligo)")
@click.option("-o",
              "--output",
              'output_file',
              required=True,
              type=click.Path(exists=False, writable=True),
              help="Pickle output file")

def generate_assignment(input_file, output_file):
    assignment = {}
    with gzip.open(input_file, "rt") as file:
        tsvreader = csv.reader(file, delimiter="\t")
        for line in tsvreader:
            oligo = line[1]
            barcode = line[0]
            barcodes_of_oligo = assignment.get(oligo,set())
            barcodes_of_oligo.add(barcode)
            assignment[oligo] = barcodes_of_oligo


    pickle.dump(assignment, open(output_file, "wb" ) )

generate_assignment()
