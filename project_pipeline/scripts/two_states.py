'''
A script to extract the proteins with two states.
'''

import os
import pandas as pd

df = pd.read_csv(snakemake.input[0], sep='\t').astype('object')

df['distinct_count'] = df.groupby('uniprot')['state'].transform('nunique')

# Dataframe with only proteins with both states
both_states = df[df['distinct_count'] == 2]

# Select only the columns we need
both_states = both_states[['uniprot', 'region_1', 'region_2']].drop_duplicates().reset_index(drop=True)

# Read which colabfold files we have
path = snakemake.input[1]
files = os.listdir(path)

file_dict = {'uniprot': [], 'filename': [], 'region_1': [], 'region_2': []}
# Create a new dataframe with the file names
for f in files:
    file_dict = {}
    uniprot = f.split('_')[0]
    region1 = df[df['uniprot'] == uniprot]['region_1']
    region2 = df[df['uniprot'] == uniprot]['region_2']

    # Populate the dictionary
    file_dict['uniprot'].append(uniprot)
    file_dict['filename'].append(f)
    file_dict['region_1'].append(region1)
    file_dict['region_2'].append(region2)

# Create a new dataframe from the dictionary
file_df = pd.DataFrame(file_dict)

# Save file
file_df.to_csv(snakemake.output[0], sep='\t', index=False)