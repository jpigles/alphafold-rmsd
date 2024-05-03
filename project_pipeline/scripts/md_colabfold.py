'''
Filter for proteins for which we have Colabfold files.
'''

import os
import pandas as pd

df = pd.read_csv(snakemake.input[0], sep='\t').astype('object')

# Read which colabfold files we have
path = snakemake.input[1]
# We have a list of directories, each with the same name as the uniprot
dirs = os.listdir(path)

file_dict = {'uniprot': [], 'filename': [], 'region_1': [], 'region_2': []}
# Create a new dataframe with the file names
for d in dirs:
    files = os.listdir(path + d)
    for f in files:
        uniprot = f.split('_')[0] # filename example: P28482_U10-000_unrelaxed_rank_001_alphafold2_multimer_v2_model_1_seed_000.pdb
        model = f.split('_')[9]
        if model != '1': # We want only the first model
            continue
        region1 = df.loc[df['uniprot'] == uniprot, 'region_1'].iloc[0]
        region2 = df.loc[df['uniprot'] == uniprot, 'region_2'].iloc[0]

        # Populate the dictionary
        file_dict['uniprot'].append(uniprot)
        file_dict['filename'].append(f)
        file_dict['region_1'].append(region1)
        file_dict['region_2'].append(region2)

# Create a new dataframe from the dictionary
file_df = pd.DataFrame(file_dict)

# Save file
file_df.to_csv(snakemake.output[0], sep='\t', index=False)