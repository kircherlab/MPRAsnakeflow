import pandas as pd
import os
import click

@click.command()
@click.option(
    "--counts",
    "counts_files",
    required=True,
    multiple=True,
    type=(str,click.Path(exists=True, readable=True)),
    help="Replicate name and assigned barcode count file",
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
    "--output",
    "output_threshold_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file.",
)
@click.option(
    "--output-threshold",
    "output_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file.",
)
def cli(counts_files, bc_thresh, output_threshold_file, output_file):
    """
    Merge the associated barcode count files of all replicates.
    """

    all_reps = []
    replicates = []
    for rep, file in counts_files:
        df = pd.read_csv(file, sep="\t")
        df['replicate'] = rep
        all_reps.append(df)
        replicates.append(rep)

    df = pd.concat(all_reps)
    df = df[df["oligo_name"] != "no_BC"]

    def pivot_table(df, replicates):
        # pivot table to make a dna and rna count column for every replicate
        df = df.pivot_table(
            values=["dna_count", "rna_count"],
            index=["barcode", "oligo_name"],
            columns="replicate",
            aggfunc='first'
        )
        df = df.sort_values("oligo_name")


        # order columns to have dna then rna count of each replicate
        col_order = sum(
            [
                ["dna_count_" + rep, "rna_count_" + rep]
                for rep in replicates
            ],
            [],
        )

        df = df.reset_index()

        df.columns = ['_'.join(col).strip() if col[1] else col[0] for col in df.columns.values]
            
        df = df[["barcode", "oligo_name"] + col_order]

        for col in col_order:
            df[col] = df[col].astype('Int32')
        
        return df

    # only keep oligo's with a number of barcodes of at least the given threshold
    df_filtered = df.groupby(["oligo_name", "replicate"]).filter(lambda x: len(x) >= bc_thresh)
    
    # write to output file
    pivot_table(df_filtered,replicates).to_csv(output_threshold_file, sep="\t", index=False, compression="gzip")
    pivot_table(df,replicates).to_csv(output_file, sep="\t", index=False, compression="gzip")


if __name__ == "__main__":
    cli()
