"""
This script checks for collisions in a fasta design file. It looks for sequences that are repeated (forward and antisense) in the design file.
Also checks for duplicated headers in the design file.

:Author: Max Schubach
:Contact: max.schubach@bih-charite.de
:Date: *19.07.2024
"""

import pyfastx
import click
import numpy as np



# options
@click.command()
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Fasta design file.",
)
@click.option(
    "--start",
    "start",
    required=True,
    type=int,
    help="Start of the sequence to look at (1 based).",
)
@click.option(
    "--length",
    "length",
    required=True,
    type=int,
    help="length of the sequence to lok at,",
)
def cli(input_file, start, length):

    seq_dict = dict()

    forward_collitions = []
    antisense_collitions = []

    # read fasta file
    fa = pyfastx.Fasta(input_file)

    # check duplicated headers
    click.echo("Searching for duplicated headers...")
    ids = fa.keys()
    if len(ids) != len(set(ids)):
        click.echo("Duplicated headers found in the design file.", err=True)
        exit(1)

    # check for illegal characters
    click.echo("Searching for illegal characters in header...")
    illegal_characters =  np.array([True if "[" in i or "]" in i else False for i in set(ids)], dtype=bool).sum()
    if illegal_characters > 0:
        click.echo(f"{illegal_characters} headers contain illegal characters ('[',']').", err=True)
        exit(1)

    # build seq dict
    click.echo("Building sequence dictionary...")
    for i in range(len(fa)):
        names = seq_dict.get(fa[i].seq, set())
        names.add(fa[i].name)
        seq_dict[fa[i].seq] = names

    # search for collisions
    click.echo("Searching for collisions...")
    for i in range(len(fa)):
        sub_seq = fa[i][start - 1 : start + length - 1]
        sub_seq_forward = sub_seq.seq
        sub_seq_antisense = sub_seq.antisense

        forward_collition = set()
        antisense_collition = set()
        for seq, names in seq_dict.items():
            if sub_seq_forward in seq:
                forward_collition.update(names)
            if sub_seq_antisense in seq:
                antisense_collition.update(names)

        if len(forward_collition) > 1:
            forward_collitions.append(forward_collition)
        if len(antisense_collition) > 1:
            antisense_collitions.append(antisense_collition)

    # unique names
    forward_collitions = [list(i) for i in set(tuple(i) for i in forward_collitions)]
    antisense_collition = [list(i) for i in set(tuple(i) for i in antisense_collition)]

    if (len(forward_collitions) > 0) or (len(antisense_collitions) > 0):
        click.echo("Collisions found:", err=True)
        for i in range(len(forward_collitions)):
            click.echo(
                f"Forward collision {i+1} for sequences:\t{"\t".join(forward_collitions[i])}", err=True
            )
        for i in range(len(antisense_collitions)):
            click.echo(
                f"Antisense collision {i+1} for sequences:\t{"\t".join(antisense_collitions[i])}",
                err=True,
            )
        exit(1)

    else:
        click.echo("No collisions found. Design file seems to be in a good shape.")


if __name__ == "__main__":
    cli()
