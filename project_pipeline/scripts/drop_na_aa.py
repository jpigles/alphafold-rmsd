import pandas as pd

df_prot = pd.read_csv(snakemake.input[0], sep='\t').astype('object')

df_prot_no_na = df_prot.dropna(axis=0, subset='PDB')

df_prot_no_na.to_csv(snakemake.output[0], sep='\t', index=False)