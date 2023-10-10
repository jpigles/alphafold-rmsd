'''
This script calculates the mean predicted aligned error for a set of 100 randomly chosen proteins
from the PDB's single-domain protein database (excluding all intrinsically disordered proteins).
See autoinhibition_protein_data/single_domain.ipynb for how these proteins were retrieved.
'''

import pandas as pd
import main

# Read in the file with the 100 randomly chosen proteins
df = pd.read_csv('./data/single_domain_uniprots.csv', sep=',').astype('object')
# Define file path for Alphafold predicted aligned error files
af_path = './data/input/Alphafold_single_domain/'

# Calculate mean PAE for each protein
df_mean_pae = main.mean_pae_single_domain(df, af_path)

# Save file
df_mean_pae.to_csv('./data/single_domain_pae.tsv', sep='\t', index=False)
