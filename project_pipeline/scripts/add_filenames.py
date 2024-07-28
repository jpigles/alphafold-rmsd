'''
A script to extract the proteins with two states.
'''

import os
import utils
import pandas as pd

'''
First we add AlphaFold2 filenames for our autoinhibited AlphaFold2 files.
'''
autoinhibited_af_path = snakemake.input[5]

df = pd.read_csv(snakemake.input[0], sep='\t').astype('object')

# Add the filenames
df_af_fn = utils.add_AF_filename(df, autoinhibited_af_path)

# Save file
df_af_fn.to_csv(snakemake.output[0], sep='\t', index=False)

'''
Then we take the proteins with two states and add the ColabFold filenames.
'''
autoinhibited_cf_path = snakemake.input[6]

df['distinct_count'] = df.groupby('uniprot')['state'].transform('nunique')

# Dataframe with only proteins with both states
both_states = df[df['distinct_count'] == 2].reset_index(drop=True)
# Drop the distinct count column
both_states = both_states.drop(columns=['distinct_count'])

# Assign state & conformation to the AlphaFold2 structure
# Find the lowest 2_comp per protein
lowest_2_comp = both_states.groupby('uniprot')['2_comp'].min().reset_index()

# Get the state, conformation, and regions for the lowest 2_comp
lowest_2_comp = pd.merge(lowest_2_comp, both_states, on=['uniprot', '2_comp'])

# To make our af to cf comparison file, keep only the uniprot, region_1, region_2, conformation, and state
af = lowest_2_comp[['uniprot', 'region_1', 'region_2', 'conformation', 'state']].drop_duplicates().reset_index(drop=True)

# Read which colabfold files we have

cf_af_df = utils.add_CF_filename(af, autoinhibited_cf_path)

# Add AlphaFold2 filenames
df_af_fn = df_af_fn[['uniprot', 'af_filename']]

cf_af_df = pd.merge(cf_af_df, df_af_fn, on='uniprot', how='left').reset_index()

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
Next, we add the AlphaFold2 filenames for our multi-domain proteins.
'''

md_af_path = snakemake.input[7]

multi = pd.read_csv(snakemake.input[1], sep='\t').astype('object')

multi_af = utils.add_AF_filename(multi, md_af_path)

# Save file
multi_af.to_csv(snakemake.output[3], sep='\t', index=False)

'''
Then we add the AlphaFold2 filenames for our obligate multi-domain proteins.
'''
obli = pd.read_csv(snakemake.input[4], sep='\t').astype('object')

obli_af = utils.add_AF_filename(obli, md_af_path)

# Save file
obli_af.to_csv(snakemake.output[6], sep='\t', index=False)

'''
Then we add filenames to our single-domain proteins.
'''
single_path = snakemake.input[8]
single = pd.read_csv(snakemake.input[3], sep='\t').astype('object')

single = single.dropna(subset=['mean_pae'])
single = single[['uniprot', 'region']]
single = utils.add_AF_filename(single, single_path)

# Save file
single.to_csv(snakemake.output[4], sep='\t', index=False)

'''
Next, we get the Colabfold filenames for our multi-domain proteins.
'''

# Read which colabfold files we have
md_cf_path = snakemake.input[9]

multi_cf_df = utils.add_CF_filename(multi, md_cf_path)

# Add AlphaFold2 filenames
multi_af = multi_af[['uniprot', 'af_filename']].drop_duplicates()

multi_cf_df = pd.merge(multi_cf_df, multi_af, on='uniprot', how='left').reset_index(drop=True)

# Save file
multi_cf_df.to_csv(snakemake.output[5], sep='\t', index=False)

'''
Get the Colabfold filenames for our 
obligate multi-domain proteins.
'''

obli_cf_path = snakemake.input[10]

obli_cf = utils.add_CF_filename(obli, obli_cf_path)

# Add AlphaFold2 filenames
obli_af = obli_af[['uniprot', 'af_filename']].drop_duplicates()

obli_cf = pd.merge(obli_cf, obli_af, on='uniprot', how='left').reset_index(drop=True)

# Save file
obli_cf.to_csv(snakemake.output[7], sep='\t', index=False)

'''
The ColabFold filenames for our 20 exemplary low-complexity and high-complexity species.
'''

species_cf_path = snakemake.input[11]

# Get rid of the distinct_counts column
df = df.drop(columns=['distinct_counts'])

species_cf_df = utils.add_CF_filename(df, species_cf_path)

species_cf_df.to_csv(snakemake.output[8], sep='\t', index=False)

'''
Make a file for AlphaFold to ColabFold comparison of species
'''

# Assign state & conformation to the AlphaFold2 structure
# Drop clusters for now
no_clust = species_cf_df.drop(columns=['cluster', 'cf_filename'])

# Find the lowest 2_comp per protein
min_2_comp = no_clust.groupby('uniprot')['2_comp'].min().reset_index()

# Take the selected state and conformation for each protein
st_conf = min_2_comp[['uniprot', 'state', 'conformation']]

# Merge them back together with main dataframe and drop unneeded info
species_cf_df = species_cf_df[['uniprot', 'region_1', 'region_2', 'cluster', 'cf_filename']].drop_duplicates().reset_index(drop=True)

species_af_df = pd.merge(species_cf_df, st_conf, on='uniprot', how='left')

# Add AlphaFold2 filenames
species_af_df = pd.merge(species_af_df, df_af_fn, on='uniprot', how='left')

species_af_df.to_csv(snakemake.output[9], sep='\t', index=False)