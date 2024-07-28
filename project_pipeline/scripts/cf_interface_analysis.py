'''
Determine the interface residues of the ColabFold-generated structures
'''
import pandas as pd
import main

path = snakemake.input[0]

df = pd.read_csv(snakemake.input[1], sep='\t')

# Get the interfaces of the proteins
print('Finding interfaces...')
interacting = main.get_af_interfaces(df, path, cluster=True)

# Save the dataframe
interacting.to_csv(snakemake.output[0], sep='\t', index=False)