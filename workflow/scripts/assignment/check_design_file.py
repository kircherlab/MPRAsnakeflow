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
import re


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
@click.option(
    "--fast-dict/--slow-string-search",
    "fast_search",
    default=True,
    help="Using a simple dictionary to find identical sequences. This is faster but uses only the whole (or center part depending on start/length) of the design file. But in theory a substring can only be present and for more correct, but slower, search use the --slow-string-search.",
)
@click.option(
    "--sequence-check",
    "sequence_check",
    type=click.Choice(["skip", "sense_only", "sense_antisense"]),
    default="sense_antisense",
    help="Choose the type of sequence check. When set to skip, the script will not check for sequence collisions. This is useful when you know collisions but still want to preoceed with the design file.",
)
@click.option(
    "--attach-sequence",
    "attach_sequence",
    required=False,
    type=(click.STRING, click.STRING),
    help="Attach a sequence left and right to each entry of the fasta file.",
)
@click.option(
    "--output",
    "output",
    required=True,
    type=click.Path(writable=True),
    help="Output reference fasqt file with attached sequences (if set). Otherwise just a copy.",
)
def cli(input_file, start, length, fast_search, sequence_check, attach_sequence, output):

    seq_dict = dict()

    forward_collitions = []
    antisense_collitions = []

    # attach sequence
    with open(output, "w") as fasta_file:
        for name, seq in pyfastx.Fasta(input_file, build_index=False):
            new_sequence = seq.upper()
            if attach_sequence:
                new_sequence = attach_sequence[0].upper() + new_sequence + attach_sequence[1].upper()
            fasta_file.write(f">{name}\n{new_sequence}\n")

    # read fasta file
    fa = pyfastx.Fasta(output)

    # check duplicated headers
    click.echo("Searching for duplicated headers...")
    ids = fa.keys()
    if len(ids) != len(set(ids)):
        click.echo("Duplicated headers found in the design file.", err=True)
        exit(1)

    # check for illegal characters
    click.echo("Searching for illegal characters in header...")
    pattern = re.compile(
        r'^[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*$'
    )
    illegal_characters_headers = [i for i in ids if pattern.fullmatch(i) is None]
    if len(illegal_characters_headers) > 0:
        click.echo(
            f"{len(illegal_characters_headers)} headers contain illegal characters. Only this expression is allowed: ^[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*$",
            err=True,
        )
        click.echo(
            f"The following headers contain illegal characters:\n" + "\n".join(illegal_characters_headers),
            err=True,
        )
        exit(1)

    # read fasta file
    if attach_sequence:
        fa = pyfastx.Fasta(output)

    if sequence_check != "skip":
        # build seq dict
        click.echo("Building sequence dictionary...")
        for i in range(len(fa)):
            if fast_search:
                seq = fa[i][start - 1 : start + length - 1].seq
            else:
                seq = fa[i].seq

            names = seq_dict.get(seq, set())
            names.add(fa[i].name)
            seq_dict[seq] = names

        # search for collisions
        click.echo("Searching for collisions...")
        for i in range(len(fa)):
            sub_seq = fa[i][start - 1 : start + length - 1]
            sub_seq_forward = sub_seq.seq
            sub_seq_antisense = sub_seq.antisense

            forward_collition = set()
            antisense_collition = set()
            if fast_search:
                forward_collition.update(seq_dict.get(sub_seq_forward, set()))
                (
                    antisense_collition.update(seq_dict.get(sub_seq_antisense, set()))
                    if sequence_check == "sense_antisense"
                    else None
                )
            else:
                for seq, names in seq_dict.items():
                    if sub_seq_forward in seq:
                        forward_collition.update(names)
                    if sequence_check == "sense_antisense" and sub_seq_antisense in seq:
                        antisense_collition.update(names)

            if len(forward_collition) > 1:
                forward_collitions.append(forward_collition)
            if len(antisense_collition) > 1:
                antisense_collitions.append(antisense_collition)

        # unique names
        forward_collitions = [list(i) for i in set(tuple(i) for i in forward_collitions)]
        antisense_collitions = [list(i) for i in set(tuple(i) for i in antisense_collitions)]

        if (len(forward_collitions) > 0) or (len(antisense_collitions) > 0):
            click.echo(
                f"Found {len(forward_collitions)} forward and {len(antisense_collitions)} antisense collisions",
                err=True,
            )
            click.echo(
                "-----------------FORWARD COLLISIONS-----------------",
                err=True,
            )
            for i in range(len(forward_collitions)):
                click.echo(
                    "\t".join(forward_collitions[i]),
                    err=True,
                )
            if sequence_check == "sense_antisense":
                click.echo(
                    "-----------------ANTISENSE COLLISIONS-----------------",
                    err=True,
                )
                for i in range(len(antisense_collitions)):
                    click.echo(
                        "\t".join(antisense_collitions[i]),
                        err=True,
                    )
            exit(1)

        else:
            click.echo("No collisions found. Design file seems to be in a good shape.")


if __name__ == "__main__":
    cli()
