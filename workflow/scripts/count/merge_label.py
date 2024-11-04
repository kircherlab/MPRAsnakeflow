# Author: Gracie Gordon 2019
# merge tables

import pandas as pd
import numpy as np
import click

# options


@click.command()
@click.option(
    "--counts",
    "counts_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Merged DNA and RNA counts.",
)
@click.option(
    "--assignment",
    "assignment_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Assignment file in tsv format.",
)
@click.option(
    "--minRNACounts",
    "minRNACounts",
    default=1,
    show_default=True,
    type=int,
    help="Minimum number of RNA counts required per barcode. If 0 pesudocounts are used.",
)
@click.option(
    "--minDNACounts",
    "minDNACounts",
    default=0,
    show_default=True,
    type=int,
    help="Minimum number of DNA counts required per barcode. If 0 pesudocounts are used.",
)
@click.option(
    "--inlclude-not-assigned-for-normalization/--exclude-not-assigned-for-normalization",
    "normalize_with_not_assigned",
    default=False,
    show_default=True,
    help="Use barcodes that are not assigned for normalization.",
)
@click.option(
    "--scaling",
    "scaling",
    default=10**6,
    show_default=True,
    type=float,
    help="Scaling parameter. Usually counts per million (10**6).",
)
@click.option(
    "--outlier-detection",
    "outlier_detection_method",
    required=False,
    type=click.Choice(['ratio_mad', 'rna_counts_zscore'], case_sensitive=False),
    help="If set it wil use obe of the outlier detection methods.",
)
@click.option(
    "--outlier-barcodes",
    "removed_barcodes_file",
    required=False,
    type=click.Path(writable=True),
    help="If outlier detection is enabled it will write out the barcodes that were removed.",
)
@click.option(
    "--outlier-ratio-mad-bins",
    "ratio_mad_n_bins",
    default = 20,
    type=int,
    help="Number of quantile bins for RNA counts where outlier removal is processed.",
)
@click.option(
    "--outlier-ratio-mad-times",
    "ratio_mad_times",
    default = 5,
    type=float,
    help="ratio diff does not be larger than the value time mad within a bin.",
)
@click.option(
    "--outlier-rna-zscore-times",
    "rna_zscore_times",
    default = 3,
    type=float,
    help="ratio diff does not be larger than the value time mad within a bin.",
)
@click.option(
    "--output",
    "output_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file.",
)
@click.option(
    "--statistic",
    "statistic_file",
    required=True,
    type=click.Path(writable=True),
    help="Statistic output file.",
)
@click.option(
    "--bcOutput",
    "bc_output_file",
    type=click.Path(writable=True),
    help="Output file for individual barcode counts",
)
def cli(
    counts_file,
    assignment_file,
    minRNACounts,
    minDNACounts,
    normalize_with_not_assigned,
    scaling,
    outlier_detection_method,
    ratio_mad_times,
    ratio_mad_n_bins,
    rna_zscore_times,
    removed_barcodes_file,
    output_file,
    statistic_file,
    bc_output_file,
):
    # pseudocount = 1 if minRNACounts == 0 or minDNACounts == 0 else 0
    pseudocountDNA = 1 if minDNACounts == 0 else 0
    pseudocountRNA = 1 if minRNACounts == 0 else 0

    # statistic
    statistic = pd.DataFrame(
        data={
            "oligos design": [0],
            "barcodes design": [0],
            "oligos dna/rna": [0],
            "barcodes dna/rna": [0],
            "matched barcodes": [0],
            "unknown barcodes dna/rna": [0],
            "% matched barcodes": [0.0],
            "total dna counts": [0],
            "total rna counts": [0],
            "avg dna counts per bc": [0.0],
            "avg rna counts per bc": [0.0],
            "barcode outlier removed": [0.0],
            "avg dna/rna barcodes per oligo": [0.0],
        }
    )

    # process association file
    click.echo("Read assignment file...")
    assoc = pd.read_csv(
        assignment_file,
        header=None,
        usecols=[0, 1],
        sep="\t",
        names=["barcode", "oligo"],
    )
    # drop duplicated barcodes!!!
    assoc.drop_duplicates("barcode", keep=False, inplace=True)
    assoc.set_index("barcode", inplace=True)

    statistic["oligos design"] = assoc.oligo.nunique()
    statistic[["barcodes design"]] = assoc.shape[0]

    # get count df
    click.echo("Read count file...")
    counts = pd.read_csv(
        counts_file, sep="\t", header=None, names=["barcode", "dna_count", "rna_count"]
    )
    # filter
    counts = counts[
        (counts["dna_count"] >= minDNACounts) & (counts["rna_count"] >= minRNACounts)
    ]

    statistic["barcodes dna/rna"] = counts.shape[0]

    # fill in labels from dictionary
    click.echo("Combine assignment with count file...")
    counts = pd.merge(assoc, counts, how="right", on="barcode")
    counts['oligo'] = counts['oligo'].fillna("no_BC")
    counts.rename(columns={"oligo": "oligo_name"}, inplace=True)
    # counts = counts[['oligo_name', 'barcode', 'dna_count', 'rna_count']]
    statistic[["unknown barcodes dna/rna"]] = counts[counts.oligo_name == "no_BC"].shape[0]

    statistic["matched barcodes"] = (
        statistic["barcodes dna/rna"] - statistic["unknown barcodes dna/rna"]
    )
    statistic["% matched barcodes"] = (
        statistic["matched barcodes"] / statistic["barcodes dna/rna"]
    ) * 100.0

    # add pseudocount if needed:
    counts["dna_count"] = counts["dna_count"] + pseudocountDNA
    counts["rna_count"] = counts["rna_count"] + pseudocountRNA


    # number of DNA and RNA counts
    total_dna_counts = sum(counts["dna_count"])
    total_rna_counts = sum(counts["rna_count"])

    unknown_dna_counts = sum(counts[counts.oligo_name == "no_BC"]["dna_count"])
    unknown_rna_counts = sum(counts[counts.oligo_name == "no_BC"]["rna_count"])

    statistic["total dna counts"] = total_dna_counts
    statistic["total rna counts"] = total_rna_counts

    statistic["avg dna counts per bc"] = (
        statistic["total dna counts"] / statistic["barcodes dna/rna"]
    )
    statistic["avg rna counts per bc"] = (
        statistic["total rna counts"] / statistic["barcodes dna/rna"]
    )

    if not normalize_with_not_assigned:
        counts = counts[counts.oligo_name != "no_BC"]
        total_dna_counts -= unknown_dna_counts
        total_rna_counts -= unknown_rna_counts

    # outlier removal

    if outlier_detection_method:
        if outlier_detection_method == "ratio_mad":
            counts, removed_barcodes = outlier_removal_by_mad(counts, ratio_mad_n_bins, ratio_mad_times)
        elif outlier_detection_method == "rna_counts_zscore":
            counts, removed_barcodes = outlier_removal_by_rna_zscore(counts, rna_zscore_times)
        else:
            raise ValueError("Outlier removal method not recognized")
    else:
        removed_barcodes = pd.DataFrame({'barcode' : []})
        
    if removed_barcodes_file:
        removed_barcodes.to_csv(removed_barcodes_file, columns=["barcode"], header=False, index=False)

    statistic["barcode outlier removed"] = removed_barcodes.shape[0]

    # BC output file
    if bc_output_file:
        counts[["barcode", "oligo_name","dna_count","rna_count"]].to_csv(bc_output_file, index=False, sep="\t", compression="gzip")

    # remove Barcorde. Not needed anymore
    counts.drop(["barcode"], axis=1, inplace=True)

    # group by oligo name 
    grouped_label = counts.groupby("oligo_name").agg(
        {"dna_count": ["sum", "count"], "rna_count": ["sum", "count"]}
    )
    
    grouped_label.reset_index(inplace=True)

    output = pd.DataFrame()

    click.echo(grouped_label.head())

    output["oligo_name"] = grouped_label["oligo_name"]
    output["dna_counts"] = grouped_label.dna_count["sum"]
    output["rna_counts"] = grouped_label.rna_count["sum"]

    statistic["oligos dna/rna"] = len(grouped_label)
    statistic["avg dna/rna barcodes per oligo"] = (
        sum(grouped_label[grouped_label["oligo_name"] != "no_BC"].dna_count["count"])
        / statistic["oligos dna/rna"]
    )


    output["dna_normalized"] = (
        ( grouped_label.dna_count["sum"]/ grouped_label.dna_count["count"] )
        / total_dna_counts
        * scaling
    )

    output["rna_normalized"] = (
        ( grouped_label.rna_count["sum"] / grouped_label.rna_count["count"] )
        / total_rna_counts
        * scaling
    )

    output["ratio"] = output["rna_normalized"] / output["dna_normalized"]
    output["log2FoldChange"] = np.log2(output.ratio)

    output["n_bc"] = grouped_label.dna_count["count"]

    click.echo(output_file)

    output.to_csv(output_file, index=False, sep="\t", compression="gzip")

    statistic.to_csv(statistic_file, index=False, sep="\t", compression="gzip")

def outlier_removal_by_rna_zscore(df, times_zscore = 3):
    df["rna_z_scores"] = df.groupby('oligo_name')['rna_count'].transform(lambda x: ( x- np.mean(x)/ np.std(x)))
    
    m = df.rna_z_scores.abs() <= times_zscore
    barcodes_removed = df[~m].barcode
    df = df[m]
    return df[m], barcodes_removed

def outlier_removal_by_mad(df, n_bins = 20, times_mad = 5):
    # df = df[df.oligo_name != "no_BC"]
    # Calculate ratio, ratio_med, ratio_diff, and mad
    df['ratio'] = np.log2(df["dna_count"] / df["rna_count"])
    df['ratio_med'] = df.groupby('oligo_name')['ratio'].transform('median')
    df['ratio_diff'] = df['ratio'] - df['ratio_med']

    # Calculate quantiles within  n_bins
    qs = np.unique(np.quantile(np.log10(df['rna_count']), np.arange(0, n_bins) / n_bins))
    
    # Create bins based on rna_count
    df['bin'] = pd.cut(np.log10(df['rna_count']), bins=qs, include_lowest=True, labels=[str(i) for i in range(0, len(qs)-1)])
    # Filter based on ratio_diff and mad
    df['mad'] = df.groupby('bin', observed=True)['ratio_diff'].transform(lambda x: np.median(np.abs(x - np.median(x)))) 
    
    m = df.ratio_diff <= times_mad * df.mad
    barcodes_removed = df[~m].barcode
    df = df[m]

    return df, barcodes_removed

if __name__ == "__main__":
    cli()
