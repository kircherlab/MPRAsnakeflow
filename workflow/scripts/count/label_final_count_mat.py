# this script processes the RNA and DNA counts and assigns the enhancer tag
#outputs dataframe

#CMD: python label_count_mat.py test.merged.H2.tsv ../lib_assoc_scripts/mp_assoc_original/bc_info_mp/Gracie_mp_filtered_coords_to_barcodes.pickle test.log2.fold.txt 

import pickle
import click
import numpy
import pandas
import sys

from Bio import SeqIO

# options
@click.command()
@click.option('--counts',
              'counts_file',
              required=True, 
              type=click.Path(exists=True, readable=True),
              help='Merged DNA and RNA counts.')
@click.option('--assignment',
              'assignment_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Assignment pickle file.')
@click.option('--design', 
              'design_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Design file in fasta format.')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file.')
def cli(counts_file, assignment_file, design_file, output_file):
    
    # process fastq
    design=open(design_file)
    fasta_dict = {rec.id : rec.seq for rec in SeqIO.parse(design, "fasta")}

    assoc=pickle.load( open(assignment_file,'rb'))

    #flip dict so that BCs are the Keys
    BC_key = {}
    for k,v in assoc.items():
        for x in v:
             BC_key.setdefault(x,k)
    
    counts=pandas.read_csv(counts_file,header='infer',sep=',')

    #fill in labels from dictionary
    label=[]
    for i in counts.Barcode:
       try:
               label.append(BC_key[i])
       except:
               label.append('no_BC')

    seqs=[]
    for l in label:
        try:
            seqs.append(str(fasta_dict[l]).upper())
        except:
            seqs.append('NA')
    counts.insert(0,'Sequence',seqs)
    counts.insert(0, 'label', label)
    click.echo(counts.head())

    counts.to_csv(output_file,header=True,index=False,sep="\t") 

if __name__ == '__main__':
    cli()

