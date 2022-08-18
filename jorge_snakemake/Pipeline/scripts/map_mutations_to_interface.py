# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 14:19:23 2020

@author: Jorge Holguin
"""

import pandas as pd
import numpy as np
import mutation_enrichment as mut
from statsmodels.stats.multitest import multipletests

def correct_multiple_testing(p_vals):
    
    for value in range(len(p_vals)):
        if p_vals[value] == 0.0:
            p_vals[value] = 1e-323
            
    # Use the Benjamini/Hochberg FDR correction algorithm
    corr_p_values = multipletests(p_vals, alpha = 0.05, method = 'fdr_bh', is_sorted=False, returnsorted=False)
    
    return corr_p_values

#=============================================================================
# Open the file as pandas dataframes
df_chunk = pd.read_csv(snakemake.input[1], sep = '\t', chunksize = 1000000, header = 0,
                      low_memory=False)

# df_chunk = pd.read_csv('../data/mutations/cosmic.tsv', sep = '\t', chunksize = 1000000, header = 0,
#                       low_memory=False)

chunk_list = []  # append each chunk df here 

cols = ['GENE_NAME', 'ACCESSION_NUMBER', 'Mutation AA', 'Mutation Description AA', 'AA_MUT_START', 'SHARED_AA', 'MUTATION_SIGNIFICANCE_TIER']

# Each chunk is in df format
for chunk in df_chunk:
    
    chunk = chunk[cols]
    chunk = chunk.loc[(chunk['Mutation Description AA'] == 'Substitution - Missense') & (chunk['MUTATION_SIGNIFICANCE_TIER'].isin(['1','2','3', 'Other']))]

    # Once the data filtering is done, append the chunk to list
    chunk_list.append(chunk)

# concat the list into dataframe, drop the unnamed column
df_cosmic = pd.concat(chunk_list)

# Re format the accession number to remove the version of the accession
df_cosmic['ACCESSION_NUMBER'] = df_cosmic['ACCESSION_NUMBER'].apply(lambda x: x.split(sep = '.')[0])

# Subset the potential drivers and potential passenger mutations
df_cosmic_patho = df_cosmic.loc[df_cosmic['MUTATION_SIGNIFICANCE_TIER'].isin(['1', '2', '3'])]
df_cosmic_other = df_cosmic.loc[df_cosmic['MUTATION_SIGNIFICANCE_TIER'].isin(['Other'])]

# Read the file with proteins and interface residues
df_prot = pd.read_csv(snakemake.input[0], sep = '\t').astype('object')

# df_prot = pd.read_csv('../data/proteins_interface.tsv', sep = '\t').astype('object')

# Keep only the rows which have a region_2 and a PDB file
df_prot = df_prot.dropna(subset = ['region_2', 'PDB ID']).reset_index(drop = True)

# Convert the interface residues from string to list
df_prot['Interface Residues'] = df_prot['Interface Residues'].apply(lambda x: [int(i) for i in x.split(sep = ',')])

# Make two copies of df_prot to store the potential drivers and the potential passenger mutations
df_prot_patho = df_prot.copy()
df_prot_other = df_prot.copy()

# Find and map potential driver mutations
df_prot_patho['Mutation_list'] = np.empty((len(df_prot_patho), 0)).tolist()
df_prot_patho['Mutation_list'] = df_prot_patho.apply(lambda x: mut.find_mutations(x['Gene_name'], x['Mutation_list'], df_cosmic_patho, 
                                                                                  identifier='GENE_NAME', col_name_mut = 'AA_MUT_START', col_name_rec='SHARED_AA'), axis = 1)

norep_patho_map_df = df_prot_patho.apply(lambda x: mut.map_mutations(x['Mutation_list'], x['Interface Residues'], x['Protein_length'], repetition = False,
                                                             region_type='interface'), axis = 1, result_type='expand')

norep_patho_map_df = norep_patho_map_df.rename({0:'mut_in_interface', 1:'mut_not_in_interface', 2:'mut_list_in_interface', 3:'mut_list_not_in_interface', 
                                    4:'interface_len', 5:'outside_len', 6:'prot_mut_rate', 7:'interface_mut_rate', 
                                    8:'outside_mut_rate', 9:'norm_interface_mut_rate', 10:'norm_outside_mut_rate', 11:'enrichment_pval'}, axis = 1)

rep_patho_map_df = df_prot_patho.apply(lambda x: mut.map_mutations(x['Mutation_list'], x['Interface Residues'], x['Protein_length'], repetition = True,
                                                             region_type='interface'), axis = 1, result_type='expand')

rep_patho_map_df = rep_patho_map_df.rename({0:'mut_in_interface', 1:'mut_not_in_interface', 2:'mut_list_in_interface', 3:'mut_list_not_in_interface', 
                                    4:'interface_len', 5:'outside_len', 6:'prot_mut_rate', 7:'interface_mut_rate', 
                                    8:'outside_mut_rate', 9:'norm_interface_mut_rate', 10:'norm_outside_mut_rate', 11:'enrichment_pval'}, axis = 1)

norep_patho_results = pd.concat([df_prot_patho, norep_patho_map_df], axis = 1).astype('object')
# Correct for multiple testing
norep_patho_results['adjusted_enrichment_pval'] = correct_multiple_testing(list(norep_patho_results['enrichment_pval']))[1]
norep_patho_results['reject_null_hypothesis'] = correct_multiple_testing(list(norep_patho_results['enrichment_pval']))[0]

rep_patho_results = pd.concat([df_prot_patho, rep_patho_map_df], axis = 1).astype('object')
# Correct for multiple testing
rep_patho_results['adjusted_enrichment_pval'] = correct_multiple_testing(list(rep_patho_results['enrichment_pval']))[1]
rep_patho_results['reject_null_hypothesis'] = correct_multiple_testing(list(rep_patho_results['enrichment_pval']))[0]

# Find and map potential passenger mutations
df_prot_other['Mutation_list'] = np.empty((len(df_prot_other), 0)).tolist()
df_prot_other['Mutation_list'] = df_prot_other.apply(lambda x: mut.find_mutations(x['Gene_name'], x['Mutation_list'], df_cosmic_other, 
                                                                                  identifier='GENE_NAME', col_name_mut = 'AA_MUT_START', col_name_rec='SHARED_AA'), axis = 1)

norep_other_map_df = df_prot_other.apply(lambda x: mut.map_mutations(x['Mutation_list'], x['Interface Residues'], x['Protein_length'], repetition = False,
                                                             region_type='interface'), axis = 1, result_type='expand')

norep_other_map_df = norep_other_map_df.rename({0:'mut_in_interface', 1:'mut_not_in_interface', 2:'mut_list_in_interface', 3:'mut_list_not_in_interface', 
                                    4:'interface_len', 5:'outside_len', 6:'prot_mut_rate', 7:'interface_mut_rate', 
                                    8:'outside_mut_rate', 9:'norm_interface_mut_rate', 10:'norm_outside_mut_rate', 11:'enrichment_pval'}, axis = 1)

rep_other_map_df = df_prot_other.apply(lambda x: mut.map_mutations(x['Mutation_list'], x['Interface Residues'], x['Protein_length'], repetition = True,
                                                             region_type='interface'), axis = 1, result_type='expand')

rep_other_map_df = rep_other_map_df.rename({0:'mut_in_interface', 1:'mut_not_in_interface', 2:'mut_list_in_interface', 3:'mut_list_not_in_interface', 
                                    4:'interface_len', 5:'outside_len', 6:'prot_mut_rate', 7:'interface_mut_rate', 
                                    8:'outside_mut_rate', 9:'norm_interface_mut_rate', 10:'norm_outside_mut_rate', 11:'enrichment_pval'}, axis = 1)

norep_other_results = pd.concat([df_prot_other, norep_other_map_df], axis = 1).astype('object')
# Correct for multiple testing
norep_other_results['adjusted_enrichment_pval'] = correct_multiple_testing(list(norep_other_results['enrichment_pval']))[1]
norep_other_results['reject_null_hypothesis'] = correct_multiple_testing(list(norep_other_results['enrichment_pval']))[0]

rep_other_results = pd.concat([df_prot_other, rep_other_map_df], axis = 1).astype('object')
# Correct for multiple testing
rep_other_results['adjusted_enrichment_pval'] = correct_multiple_testing(list(rep_other_results['enrichment_pval']))[1]
rep_other_results['reject_null_hypothesis'] = correct_multiple_testing(list(rep_other_results['enrichment_pval']))[0]

# Drop some columns that will not be used in later steps
norep_patho_results = norep_patho_results.drop(['prot_mut_rate', 'interface_mut_rate', 'outside_mut_rate', 
                                                'norm_interface_mut_rate', 'norm_outside_mut_rate'], axis = 'columns')

rep_patho_results = rep_patho_results.drop(['prot_mut_rate', 'interface_mut_rate', 'outside_mut_rate', 
                                                'norm_interface_mut_rate', 'norm_outside_mut_rate'], axis = 'columns')

norep_other_results = norep_other_results.drop(['prot_mut_rate', 'interface_mut_rate', 'outside_mut_rate', 
                                                'norm_interface_mut_rate', 'norm_outside_mut_rate'], axis = 'columns')

rep_other_results = rep_other_results.drop(['prot_mut_rate', 'interface_mut_rate', 'outside_mut_rate', 
                                                'norm_interface_mut_rate', 'norm_outside_mut_rate'], axis = 'columns')

# Convert list to strings 
def to_string(x):
    x = str(x)
    
    if x == 'nan':
        return np.nan
    elif '{' in x:
        x = x.replace('{', '').replace('}', '').replace(' ', '')
        return x
    elif '[' in x:
        x = x.replace('[', '').replace(']', '').replace(' ', '')
        return x

cols = ['Interface Residues', 'Mutation_list', 'mut_list_in_interface', 'mut_list_not_in_interface']

for col in cols:
    
    norep_patho_results.loc[:, col] = norep_patho_results[col].apply(to_string)
    rep_patho_results.loc[:, col] = rep_patho_results[col].apply(to_string)
    
    norep_other_results.loc[:, col] = norep_other_results[col].apply(to_string)
    rep_other_results.loc[:, col] = rep_other_results[col].apply(to_string)

# Save the files
norep_patho_results.to_csv(snakemake.output[0], sep = '\t', index = False)
rep_patho_results.to_csv(snakemake.output[1], sep = '\t', index = False)

# norep_other_results.to_csv(snakemake.output[2], sep = '\t', index = False)
# rep_other_results.to_csv(snakemake.output[3], sep = '\t', index = False)



























