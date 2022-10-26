import click
from common import read_fastq
import gzip


def read_sequence_files(read_file, bc_file):
    for read, bc in zip(read_fastq(read_file), read_fastq(bc_file)):
        seqid_read, seq_read, qual_read = read
        seqid_read = seqid_read.split(" ")[0]
        seqid_bc, seq_bc, qual_bc = bc
        seqid_bc = seqid_bc.split(" ")[0]
        if seqid_read != seqid_bc:
            raise Exception('Sequence IDs do not match: %s != %s' % (seqid_read, seqid_bc))
        seqid = "%s XI:Z:%s,YI:Z:%s" % (seqid_read, seq_bc, qual_bc)
        yield seqid, seq_read, qual_read
    return


@click.command()
@click.option('--reads', '-r',
              "read_file",
              type=click.Path(exists=True, readable=True),
              required=True)
@click.option('--barcodes', '-b',
              "barcode_file",
              type=click.Path(exists=True, readable=True),
              required=True)
def cli(read_file, barcode_file):
    
    with gzip.open(read_file, 'rt') as r_file, gzip.open(barcode_file, 'rt') as bc_file:
        for seqid, seq, qual in read_sequence_files(r_file, bc_file):
            click.echo("@%s\n%s\n+\n%s" % (seqid, seq, qual))


if __name__ == '__main__':
    cli()
