'''
This script calculates the mean predicted aligned error for the subset of proteins labeled as "single domain" by
SCOP-e on the PDB that have only one domain annotated on UniProt and the subset of proteins labeled as "multi domain" by
SCOP-e on the PDB that have multiple domains annotated on Uniprot. 
See autoinhibition_protein_data/single_domain.ipynb for how these proteins were retrieved.
'''

import pandas as pd
import main

# Read in the file for the single domain proteins
df = pd.read_csv('./data/single_domain_domains.csv', sep=',').astype('object')
# Define file path for Alphafold predicted aligned error files
af_path = './data/input/Alphafold_single_domain/'

# Calculate mean PAE for each protein
df_mean_pae = main.mean_pae_single_domain(df, af_path)

# Save file
df_mean_pae.to_csv('./data/single_domain_pae.tsv', sep='\t', index=False)

# Read in the file for the multi domain proteins
df = pd.read_csv('./data/multi_domain_domains.csv', sep=',').astype('object')

# Define file path for Alphafold predicted aligned error files
af_path = './data/input/Alphafold_multi_domain/'

# Calculate mean PAE for each protein
df_mean_pae = main.mean_pae_single_domain(df, af_path)

# Save file
df_mean_pae.to_csv('./data/multi_domain_pae.tsv', sep='\t', index=False)
