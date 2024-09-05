import pandas as pd
import os
import click

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
		df = pd.read_csv(file, sep="\t")
		df['replicate'] = curr_rep
		all_reps.append(df)

	df = pd.concat(all_reps)
	df = df[df["name"] != "no_BC"]

	# only keep oligo's with a number of barcodes of at least the given threshold
	df_filtered = df.groupby(["name", "replicate"]).filter(lambda x: len(x) >= bc_thresh)

	# pivot table to make a dna and rna count column for every replicate
	df_filtered = df_filtered.pivot_table(
		values=["dna_count", "rna_count"],
		index=["Barcode", "name"],
		columns="replicate",
		aggfunc='first'
	)
	df_filtered = df_filtered.sort_values("name")


	# order columns to have dna then rna count of each replicate
	col_order = sum(
		[
			["dna_count_" + rep, "rna_count_" + rep]
			for rep in replicates
		],
		[],
	)

	df_filtered = df_filtered.reset_index()

	df_filtered.columns = ['_'.join(col).strip() if col[1] else col[0] for col in df_filtered.columns.values]
		
	df_filtered = df_filtered[["Barcode", "name"] + col_order]

	for col in col_order:
		df_filtered[col] = df_filtered[col].astype('Int64')

	# write to output file
	df_filtered.to_csv(output_file, sep="\t", index=False, compression="gzip")


if __name__ == "__main__":
	cli()
