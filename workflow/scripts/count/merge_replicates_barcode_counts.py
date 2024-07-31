import polars as pl
import os
import click
import gzip


@click.command()
@click.option(
    "--counts",
    "counts_files",
    required=True,
    multiple=True,
    type=click.Path(exists=True, readable=True),
    help="Assigned barcode count file",
)
@click.option(
    "--threshold",
    "bc_thresh",
    required=False,
    default=10,
    type=int,
    help="Number of required barcodes (default 10)",
)
@click.option(
    "--replicate",
    "replicates",
    multiple=True,
    type=str,
    help="replicate name",
    required=True,
)
@click.option(
    "--output",
    "output_file",
    required=True,
    default="/data/humangen_kircherlab/pia/mprasnakeflow_test/count_basic/test_output.tsv.gz",
    type=click.Path(writable=True),
    help="Output file.",
)
def cli(counts_files, bc_thresh, replicates, output_file):
    """
    Merge the associated barcode count files of all replicates.
    """

    # ensure there are as many replicates as there are files
    if len(replicates) != len(counts_files):
        raise (
            click.BadParameter(
                "Number of replicates ({}) doesn't equal the number of files ({}).".format(
                    len(replicates), len(counts_files)
                )
            )
        )

    # check if every file exists
    for file in counts_files:
        if not os.path.exists(file):
            raise (click.BadParameter("{}: file not found".format(file)))

    all_reps = []
    for file in counts_files:
        curr_rep = -1
        # find the replicate name of the current file
        for rep in replicates:
            if rep in os.path.basename(file).split("_")[1]:
                curr_rep = rep
                break
        if curr_rep == -1:
            raise (click.BadParameter("{}: incorrect file".format(file)))
        all_reps.append(
            pl.read_csv(file, separator="\t").with_columns(
                pl.lit(curr_rep).alias("replicate")
            )
        )

    df = pl.concat(all_reps)
    df = df.filter(pl.col("name") != "no_BC")
    # only keep oligo's with a number of barcodes of at least the given threshold
    df_filtered = df.filter(
        pl.col("Barcode").count().over(["name", "replicate"]) >= bc_thresh
    )

    # pivot table to make a dna and rna count column for every replicate
    df_filtered = df_filtered.pivot(
        values=["dna_count", "rna_count"],
        index=["Barcode", "name"],
        on="replicate",
    )
    df_filtered = df_filtered.sort("name")

    # order columns to have dna then rna count of each replicate
    col_order = sum(
        [
            ["dna_count_" + rep, "rna_count_" + rep]
            for rep in replicates
        ],
        [],
    )
    df_filtered = df_filtered.select(["Barcode", "name"] + col_order)

    with gzip.open(output_file, "wb") as f:
        df_filtered.write_csv(f, separator="\t")


if __name__ == "__main__":
    cli()
