import click
import pandas as pd
import numpy as np
from mpralib.mpradata import MPRABarcodeData
from Bio import SeqIO
import json


@click.group(help="Command line interface to generate quality metrics from MPRAsnakeflow.")
def cli():
    pass


@cli.command(help="QC metrics of MPRAsnakeflow experiment workflow")
@click.option(
    "--assignment", "assignment_file", type=click.Path(exists=True), required=True, help="Path to the assignment file."
)
@click.option("--design", "design_file", type=click.Path(exists=True), required=True, help="Path to the design file.")
@click.option("--output", "output_file", type=click.Path(), help="Path to the output file for the metrics (json file).")
def assignment(assignment_file, design_file, output_file):

    output = {}
    # median_assigned_barocdes
    assignment_df = pd.read_csv(assignment_file, sep="\t", header=None)

    assignment_grouped = assignment_df.groupby(1).size()
    output["median_assigned_barcodes"] = int(assignment_grouped.median())

    # fraction_assigned_oligos
    assigned_oligos = len(assignment_grouped)
    total_oligos = count_fasta_entries(design_file)

    output["fraction_assigned_oligos"] = round(assigned_oligos / total_oligos, 4)

    # output
    json_string = json.dumps(output, indent=4)
    click.echo(json_string)

    if output_file:
        with open(output_file, "w") as f:
            f.write(json_string)


@cli.command(help="QC metrics of MPRAsnakeflow experiment workflow.")
@click.option("--barcode", "barcode_file", type=click.Path(exists=True), required=True, help="Path to the barcode count file.")
@click.option(
    "--assignment", "assignment_file", type=click.Path(exists=True), required=True, help="Path to the assignment file."
)
@click.option(
    "--bc-threshold",
    "bc_threshold",
    required=False,
    default=10,
    type=int,
    help="Using a barcode threshold for output correlations.",
)
@click.option("--output", "output_file", type=click.Path(), help="Path to the output file for the metrics.")
def experiment(barcode_file, assignment_file, bc_threshold, output_file):

    output = {}

    mpra_oligo_data = MPRABarcodeData.from_file(barcode_file).oligo_data

    # median_barcodes_passing_filtering
    output["median_barcodes_passing_filtering"] = int(median_barcodes_passing_filtering(mpra_oligo_data))

    # median_rna_read_count
    output["median_rna_read_count"] = int(median_rna_read_count(mpra_oligo_data))

    # now with threshold
    mpra_oligo_data.barcode_threshold = bc_threshold

    # pearson_correlation
    output["pearson_correlation"] = float(pearson_correlation(mpra_oligo_data).round(4))

    # fraction_oligos_passing
    assignment_df = pd.read_csv(assignment_file, sep="\t", header=None)
    assignment_grouped = assignment_df.groupby(1).size()
    assigned_oligos = len(assignment_grouped)
    output["fraction_oligos_passing"] = float(round(fraction_oligos_passing(mpra_oligo_data, assigned_oligos), 4))

    # output
    json_string = json.dumps(output, indent=4)
    click.echo(json_string)

    if output_file:
        with open(output_file, "w") as f:
            f.write(json_string)


def fraction_oligos_passing(mpra_oligo_data, assigned_oligos) -> float:
    n_oligos_replicate = []
    for replicate in mpra_oligo_data.obs_names:
        replicate_data = mpra_oligo_data.data[replicate, :]
        replicate_data = replicate_data[:, replicate_data.layers["barcode_counts"] >= mpra_oligo_data.barcode_threshold]
        n_oligos_replicate += [len(replicate_data.var["oligo"])]

    return np.median(n_oligos_replicate) / assigned_oligos


def pearson_correlation(mpra_oligo_data) -> float:

    return np.median(np.triu(mpra_oligo_data.correlation()))


def count_fasta_entries(fasta_file) -> int:
    count = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        count += 1
    return count


def median_barcodes_passing_filtering(mpra_oligo_data) -> int:

    n_barcodes_replicate = []
    for replicate in mpra_oligo_data.obs_names:
        replicate_data = mpra_oligo_data.data[replicate, :]
        replicate_data = replicate_data[:, replicate_data.layers["barcode_counts"] >= mpra_oligo_data.barcode_threshold]
        n_barcodes_replicate += [np.median(replicate_data.layers["barcode_counts"])]

    return int(np.median(n_barcodes_replicate))


def median_rna_read_count(mpra_oligo_data) -> int:
    n_rna_replicate = []
    for replicate in mpra_oligo_data.obs_names:
        replicate_data = mpra_oligo_data.data[replicate, :]
        replicate_data = replicate_data[:, replicate_data.layers["barcode_counts"] >= mpra_oligo_data.barcode_threshold]
        n_rna_replicate += [np.median(replicate_data.layers["rna"])]

    return int(np.median(n_rna_replicate))


if __name__ == "__main__":
    cli()
