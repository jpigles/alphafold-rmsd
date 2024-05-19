'''
A script to extract the proteins with two states.
'''

import os
import utils
import pandas as pd

'''
First we add AlphaFold2 filenames for our autoinhibited AlphaFold2 files.
'''
autoinhibited_af_path = snakemake.input[4]

df = pd.read_csv(snakemake.input[0], sep='\t').astype('object')

# Add the filenames
df_af_fn = utils.add_AF_filename(df, autoinhibited_af_path)

# Save file
df_af_fn.to_csv(snakemake.output[0], sep='\t', index=False)

'''
Then we take the proteins with two states and add the ColabFold filenames.
'''
autoinhibited_cf_path = snakemake.input[5]

df['distinct_count'] = df.groupby('uniprot')['state'].transform('nunique')

# Dataframe with only proteins with both states
both_states = df[df['distinct_count'] == 2]

# Assign state & conformation to the AlphaFold2 structure
# Find the lowest 2_comp per protein
lowest_2_comp = both_states.groupby('uniprot')['2_comp'].min().reset_index()

# Get the state, conformation, and regions for the lowest 2_comp
lowest_2_comp = pd.merge(lowest_2_comp, both_states, on=['uniprot', '2_comp'])

# To make our af to cf comparison file, keep only the uniprot, region_1, region_2, conformation, and state
af = lowest_2_comp[['uniprot', 'region_1', 'region_2', 'conformation', 'state']].drop_duplicates().reset_index(drop=True)

# Read which colabfold files we have

cf_af_df = utils.add_CF_filename(af, autoinhibited_cf_path)

# Save file
cf_af_df.to_csv(snakemake.output[1], sep='\t', index=False)

'''
Then we make our autoinhibited two-state ColabFold to PDB comparison file.
'''

# Get chain and pdb information
pdb = pd.read_csv(snakemake.input[2], sep='\t').astype('object')

# For pdb, keep only uniprot, pdb, and chain
pdb = pdb[['uniprot', 'pdb', 'chain']]

# Get ColabFold filenames and clusters
df_cf_fn = utils.add_CF_filename(both_states, autoinhibited_cf_path)

# Add chain info
merged = pd.merge(df_cf_fn, pdb, on=['uniprot', 'pdb'], how='left')

# Save file
merged.to_csv(snakemake.output[2], sep='\t', index=False)

'''
Next, we add the AlphaFold2 filenames for our multi-domain proteins.'''

md_af_path = snakemake.input[6]

multi = pd.read_csv(snakemake.input[1], sep='\t').astype('object')

multi_af = utils.add_AF_filename(multi, md_af_path)

# Save file
multi_af.to_csv(snakemake.output[3], sep='\t', index=False)

'''
Then we add filenames to our single-domain proteins.
'''
single_path = snakemake.input[7]
single = pd.read_csv(snakemake.input[3], sep='\t').astype('object')

single = single.dropna(subset=['mean_pae'])
single = single[['uniprot', 'region']]
single = utils.add_AF_filename(single, single_path)

# Save file
single.to_csv(snakemake.output[4], sep='\t', index=False)

'''
Lastly, we get the Colabfold filenames for our multi-domain proteins.
'''

# Read which colabfold files we have
md_cf_path = snakemake.input[8]

file_df = utils.add_CF_filename(multi, md_cf_path)

# Save file
file_df.to_csv(snakemake.output[5], sep='\t', index=False)