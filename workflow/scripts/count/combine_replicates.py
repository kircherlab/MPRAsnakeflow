# Author max Schubach, 2023
# count tables of replicates into one table using sum of counts and also the mean ratio. Add the labelsof oligos if present.

import pandas as pd
import numpy as np
import click

# options


@click.command()
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Final master table with all replicates in it.",
)
@click.option(
    "--labels",
    "label_file",
    required=False,
    type=click.Path(exists=True, readable=True),
    help="Label tsv file with labels of oligos.",
)
@click.option(
    "--output",
    "output_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file.",
)
def cli(input_file, label_file, output_file):

    df_allreps = pd.read_csv(input_file, sep="\t")

    total_dna_counts = df_allreps[["dna_counts"]].sum().iloc[0]
    total_rna_counts = df_allreps[["rna_counts"]].sum().iloc[0]

    df_out = combine_replicates(label_file, df_allreps, total_dna_counts, total_rna_counts)

    df_out.to_csv(output_file, sep="\t", index=False)

def combine_replicates(label_file, df_allreps, total_dna_counts, total_rna_counts):
    df_allreps = df_allreps.groupby(by=["condition", "name"]).aggregate(
        {
            "replicate": "count",
            "dna_counts": ["sum", "mean"],
            "rna_counts": ["sum", "mean"],
            "dna_normalized": "mean",
            "rna_normalized": "mean",
            "ratio": "mean",
            "log2": "mean",
            "n_obs_bc": ["sum", "mean"],
        }
    )
    
    df_allreps = df_allreps.reset_index()
    df_out = df_allreps.iloc[:, 0:2]
    df_out.columns = ["condition", "name"]

    df_out["replicates"] = df_allreps.replicate["count"]

    scaling = 10**6

    df_out["dna_counts"] = df_allreps.dna_counts["sum"]
    df_out["rna_counts"] = df_allreps.rna_counts["sum"]


    df_out["dna_normalized"] = df_out["dna_counts"] / total_dna_counts * scaling
    df_out["rna_normalized"] = df_out["rna_counts"] / total_rna_counts * scaling

    df_out["ratio"] = df_out["rna_normalized"] / df_out["dna_normalized"]
    df_out["log2"] = np.log2(df_out.ratio)

    df_out["mean_dna_counts"] = df_allreps.dna_counts["mean"]
    df_out["mean_rna_counts"] = df_allreps.rna_counts["mean"]
    df_out["mean_dna_normalized"] = df_allreps.dna_normalized["mean"]
    df_out["mean_rna_normalized"] = df_allreps.rna_normalized["mean"]
    df_out["mean_ratio"] = df_allreps.ratio["mean"]
    df_out["mean_log2"] = df_allreps.log2["mean"]
    
    df_out["mean_n_obs_bc"] = df_allreps.n_obs_bc["mean"]


    if label_file:
        df_labels = pd.read_csv(label_file, sep="\t", names=["name", "Label"])
        df_out = df_out.join(df_labels.set_index("name"), on="name")
    return df_out


if __name__ == "__main__":
    cli()
