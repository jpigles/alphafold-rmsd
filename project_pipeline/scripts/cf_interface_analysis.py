'''
Determine the interface residues of the ColabFold-generated structures
'''
import pandas as pd
import main

auto_path = snakemake.input[0]
obli_path = snakemake.input[1]

auto_df = pd.read_csv(snakemake.input[2], sep='\t')
obli_df = pd.read_csv(snakemake.input[3], sep='\t')

# Get the interfaces of the autoinhibited proteins
print('Finding autoinhibited interfaces...')
auto_interacting = main.get_af_interfaces(auto_df, auto_path, cluster=True)

# Get the interfaces of the obligate proteins
print('Finding obligate interfaces...')
obli_interacting = main.get_af_interfaces(obli_df, obli_path, cluster=True)

# Save the dataframes
auto_interacting.to_csv(snakemake.output[0], sep='\t', index=False)

obli_interacting.to_csv(snakemake.output[1], sep='\t', index=False)