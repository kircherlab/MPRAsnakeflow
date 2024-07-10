import pandas as pd

df_all = pd.read_csv("assigned_counts.default.tsv.gz", sep='\t', compression='gzip')
df_filter = pd.read_csv("assigned_counts.default.tsv.gz", sep='\t', compression='gzip')

new_row = pd.DataFrame({'Unnamed: 0': ['Sequences in design file'], 'Counts': int(300)})

# Concatenate the new row DataFrame with the original DataFrame
df_all = pd.concat([new_row, df_all]).reset_index(drop=True).rename(columns={'Counts': 'Before filter', 'Unnamed: 0': ' '})
df_filter = (
    pd.concat([new_row, df_filter])
    .reset_index(drop=True)
    .rename(columns={"Counts": "Before filter", "Unnamed: 0": " "})
)
df = df_all.merge(df_filter, on=' ')

df
