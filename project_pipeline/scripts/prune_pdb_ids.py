'''
Created on Friday Oct 14th 2022

Author: Brooks Perkins-Jechow

'''

import pandas as pd

def prune_extra_chains(pdb_ids):
    
    #Create list of ids without chains.
    pdb_ids_only = []
    single_chain_ids = ''

    #Collect the IDs and add them to pdb_ids_only
    for pdb_id in pdb_ids:
        id_only = pdb_id[:4]
        lowercase_id = id_only.lower()
        pdb_ids_only.append(lowercase_id)

    #Determine if a PDB ID is unique. If it isn't unique, remove it.
    for pdb_id in pdb_ids_only:
        if pdb_ids_only.count(pdb_id) != 1:
            pdb_ids_only.remove(pdb_id)

    #If it is unique, add io 
        else:
            single_chain_ids = single_chain_ids + pdb_id + ' '
    return single_chain_ids.strip()


df_prot = pd.read_csv(snakemake.input[0], sep='\t').astype('object')

#Curate the list of PDB IDs
for i in range(len(df_prot)):
    if type(df_prot.loc[i, 'PDB']) == float:
        df_prot = df_prot.drop(index=[i])
    else:
        pdb_ids_unpruned = df_prot.loc[i, 'PDB']
        pdb_ids_list = pdb_ids_unpruned.strip().split(sep = ' ')
        pdb_ids_pruned = prune_extra_chains(pdb_ids_list)
        df_prot.loc[i, 'PDB'] = pdb_ids_pruned

df_prot.to_csv(snakemake.output[0], sep='\t', index=False)