import pandas as pd
import main

df = pd.read_csv('.data/input/rmsds.tsv', sep='\t').astype('object')

# Categorize proteins that have both open and closed structures
two_conf = main.two_state_proteins(df)

# Filter df to proteins in two_conf
df_two_conf = df[df['UniProt'].isin(two_conf)]

# Calculate the disorder for region 1 of each protein of interest.
df_disorder_1 = main.calculate_disorder(df_two_conf)

# Remove unneeded columns
df_disorder_1 = df_disorder_1[['uniprot', 'pdb', 'complex_rmsd', '2_aligned', '2_comp', 'state', 'conformation', 'percent_disorder_1']]

# Calculate percent of structures within 2.5A of closed conformation
df_disorder_1['2_comp'] = pd.to_numeric(df_disorder_1['2_comp'])
df_closed_2 = df_disorder_1[(df_disorder_1['Conformation'] == 'Closed') & (df_disorder_1['2_comp'] <= 2.5)]
percent_closed = len(df_closed_2)/len(df_disorder_1)*100
print(percent_closed)

df_open = df_disorder_1[(df_disorder_1['Conformation'] == 'Open') & (df_disorder_1['2_comp'] > 2.5)]
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