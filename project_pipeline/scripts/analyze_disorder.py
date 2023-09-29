'''
This script takes in the file with classified proteins (i.e. proteins with assigned states and conformations), determines which ones
have both conformations, and then calculates several values. The first is the disorder of region 1, which allows us to classify
that region as 'structured' or 'disordered'. We then calculate three values: the average predicted aligned error for each region
compared against itself, and then both regions compared to each other.
'''

import pandas as pd
import main

df = pd.read_csv('./data/classified_files_3.tsv', sep='\t').astype('object')
af_path = './data/input/Alphafold_cif/'

# # Categorize proteins that have both open and closed structures
# two_conf = main.two_state_proteins(df)

# # Filter df to proteins in two_conf
# df_two_conf = df[df['uniprot'].isin(two_conf)].reset_index(drop=True)

# Calculate the disorder for region 1 of each protein of interest.
df_disorder_1 = main.calculate_disorder(df)


# Calculate percent of structures within 2.5A of closed conformation
df_disorder_1['2_comp'] = pd.to_numeric(df_disorder_1['2_comp'])
df_closed_2 = df_disorder_1[(df_disorder_1['conformation'] == 'Closed') & (df_disorder_1['2_comp'] <= 2.5)]
percent_closed = len(df_closed_2)/len(df_disorder_1)*100
print(percent_closed)

df_open = df_disorder_1[(df_disorder_1['conformation'] == 'Open') & (df_disorder_1['2_comp'] > 2.5)]
percent_open = len(df_open)/len(df_disorder_1)*100
print(percent_open)

# Separate ARs that are mainly structured from disordered ones
df_disorder_1['percent_disorder_1'] = pd.to_numeric(df_disorder_1['percent_disorder_1'])
for i in range(len(df_disorder_1)):
    if df_disorder_1.loc[i, 'percent_disorder_1'] < 50:
        df_disorder_1.loc[i, 'ar_disorder'] = 'structured'
    else:
        df_disorder_1.loc[i, 'ar_disorder'] = 'disordered'

# Correlate RMSD of AR and AlphaFold confidence scores for the AR. 
# plDDT score per residue is given under _ma_qa_metric_local.metric_value
df_mean_plddt = main.mean_plddt(df_disorder_1, af_path)

# Calculate mean predicted aligned error for each region compared against itself and the other region
df_mean_pae = main.mean_paes(df_mean_plddt, af_path)

# Save file
df_mean_pae.to_csv('./data/disorder.tsv', sep='\t', index=False)