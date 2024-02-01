'''
A script to extract the proteins with two states.
'''

import os
import pandas as pd

df = pd.read_csv(snakemake.input[0], sep='\t').astype('object')

df['distinct_count'] = df.groupby('uniprot')['state'].transform('nunique')

# Dataframe with only proteins with both states
both_states = df[df['distinct_count'] == 2]

# Find the lowest 2_comp per protein
lowest_2_comp = both_states.groupby('uniprot')['2_comp'].min().reset_index()

# Get the state, conformation, and regions for the lowest 2_comp
lowest_2_comp = pd.merge(lowest_2_comp, both_states, on=['uniprot', '2_comp'])

cols = lowest_2_comp.columns.tolist()
cols = cols[:1] + cols[2:9] + cols[1:2] + cols[9:]
lowest_2_comp = lowest_2_comp[cols]

# To make our af to cf comparison file, keep only the uniprot, region_1, region_2, conformation, and state
af = lowest_2_comp[['uniprot', 'region_1', 'region_2', 'conformation', 'state']].drop_duplicates().reset_index(drop=True)

# Read which colabfold files we have
path = snakemake.input[1]
files = os.listdir(path)

file_dict = {'uniprot': [], 'cluster': [], 'filename': [], 'region_1': [], 'region_2': [], 'conformation': [], 'state': []}
# Create a new dataframe with the file names
for f in files:
    uniprot = f.split('_')[0] # filename example: P28482_U10-000_unrelaxed_rank_001_alphafold2_multimer_v2_model_1_seed_000.pdb
    cluster = f.split('_')[1]
    model = f.split('_')[9]
    if model != '1': # We want only the first model
        continue
    region1 = af.loc[af['uniprot'] == uniprot, 'region_1'].iloc[0]
    region2 = af.loc[af['uniprot'] == uniprot, 'region_2'].iloc[0]
    conformation = af.loc[af['uniprot'] == uniprot, 'conformation'].iloc[0]
    state = af.loc[af['uniprot'] == uniprot, 'state'].iloc[0]

    # Populate the dictionary
    file_dict['uniprot'].append(uniprot)
    file_dict['cluster'].append(cluster)
    file_dict['filename'].append(f)
    file_dict['region_1'].append(region1)
    file_dict['region_2'].append(region2)

# Create a new dataframe from the dictionary
cf_af_df = pd.DataFrame(file_dict)

# Save file
cf_af_df.to_csv(snakemake.output[0], sep='\t', index=False)

# Make our pdb to cf comparison file

# Get chain and pdb information
pdb = pd.read_csv(snakemake.input[2], sep='\t').astype('object')

# For pdb, keep only uniprot, pdb, and chain
pdb = pdb[['uniprot', 'pdb', 'chain']]

# Merge the dataframes so that we have rows for every cluster-pdb pair within the same uniprot
merged = pd.merge(cf_af_df, pdb, on='uniprot')

# Put pdb in the third column
cols = merged.columns.tolist()
cols = cols[:2] + cols[-2:-1] + cols[2:-2] + cols[-1:]
merged = merged[cols]

# Save file
merged.to_csv(snakemake.output[1], sep='\t', index=False)