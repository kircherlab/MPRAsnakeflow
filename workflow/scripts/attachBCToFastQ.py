import click
from common import read_fastq
import gzip
from typing import Optional


def read_sequence_files(
    read_file: str,
    bc_file: str,
    use_BC_reverse_complement: bool = False,
    add_sequence_left: Optional[str] = None,
    add_sequence_right: Optional[str] = None,
):
    for read, bc in zip(read_fastq(read_file), read_fastq(bc_file)):
        seqid_read, seq_read, qual_read = read

        if add_sequence_left:
            seq_read = add_sequence_left + str(seq_read)
            qual_read = "I" * len(add_sequence_left) + str(qual_read)
        if add_sequence_right:
            seq_read = str(seq_read) + add_sequence_right
            qual_read = str(qual_read) + "I" * len(add_sequence_right)

        seqid_read = str(seqid_read).split(" ")[0]
        seqid_bc, seq_bc, qual_bc = bc
        seqid_bc = str(seqid_bc).split(" ")[0]
        if seqid_read != seqid_bc:
            raise Exception("Sequence IDs do not match: %s != %s" % (seqid_read, seqid_bc))
        if use_BC_reverse_complement:
            seq_bc = reverse_complement(seq_bc)
            qual_bc = str(qual_bc)[::-1]
        seqid = "%s XI:Z:%s,YI:Z:%s" % (seqid_read, seq_bc, qual_bc)
        if len(seqid) >= 255:
            """255 character limit for QNAME exceeded, No BC"""
            seqid = "%s XI:Z:%s,YI:Z:%s" % (seqid_read, "N", "I")

        yield seqid, seq_read, qual_read
    return


@click.command()
@click.option("--reads", "-r", "read_file", type=click.Path(exists=True, readable=True), required=True)
@click.option("--barcodes", "-b", "barcode_file", type=click.Path(exists=True, readable=True), required=True)
@click.option("--reverse-complement", "use_reverse_complement", is_flag=True)
@click.option(
    "--attach-sequence", "attach_sequence", required=False, type=(click.Choice(["left", "right", "both"]), click.STRING)
)
def cli(read_file, barcode_file, use_reverse_complement, attach_sequence):

    with gzip.open(read_file, "rt") as r_file, gzip.open(barcode_file, "rt") as bc_file:
        inputs = {"read_file": r_file, "bc_file": bc_file, "use_BC_reverse_complement": use_reverse_complement}
        if attach_sequence:
            if attach_sequence[0] == "left" or attach_sequence[0] == "both":
                inputs["add_sequence_left"] = attach_sequence[1]
            if attach_sequence[0] == "right" or attach_sequence[0] == "both":
                inputs["add_sequence_right"] = reverse_complement(attach_sequence[1])

        for seqid, seq, qual in read_sequence_files(**inputs):
            click.echo("@%s\n%s\n+\n%s" % (seqid, seq, qual))


def reverse_complement(seq):
    complementary = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(reversed([complementary[i] for i in seq]))


if __name__ == "__main__":
    cli()
