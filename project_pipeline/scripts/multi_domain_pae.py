'''
This script calculates the same three predicted aligned error values calculated for our autoinhibited protein set but for 
multi-domain proteins (as classified by the PDB) with two regions.
'''

import pandas as pd
import main

df = pd.read_csv('./data/multi_domain_regions.csv').astype('object')
af_path = './data/input/Alphafold_multi_domain/'
affix = 'AF-'
suffix = '-F1-predicted_aligned_error_v4.json'

# Calculate mean predicted aligned error for each region compared against itself and the other region
df_mean_pae = main.mean_paes(df, af_path, affix, suffix)

# Save file
df_mean_pae.to_csv('./data/multi_domain_pae.tsv', sep='\t', index=False)