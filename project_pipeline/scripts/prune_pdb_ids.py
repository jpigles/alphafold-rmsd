'''
Created on Friday Oct 14th 2022

Author: Brooks Perkins-Jechow

'''

import pandas as pd

def prune_extra_chains(pdb_ids):
    pdb_ids_only = []
    single_chain_ids = []
    for pdb_id in pdb_ids:
        id_only = pdb_id[:4]
        lowercase_id = id_only.lower()
        pdb_ids_only.append(lowercase_id)
    for pdb_id in pdb_ids_only:
        if pdb_ids_only.count(pdb_id) != 1:
            pdb_ids_only.remove(pdb_id)
        else:
            single_chain_ids.append(pdb_id)
    return single_chain_ids


df_prot = pd.read_csv(snakemake.input[0], sep='\t').astype('object')

#Remove any empty PDB ID entries
for i in range(len(df_prot)):
    if type(df_prot.loc[i, 'PDB']) == float:
        df_prot.drop(index=i)

for i in range(len(df_prot)):
    pdb_ids_unpruned = df_prot.loc[i, 'PDB']
    pdb_ids_list = pdb_ids_unpruned.strip().split(sep = ' ')
    pdb_ids_pruned = prune_extra_chains(pdb_ids_list)
    df_prot.loc[i, 'PDB'] = str(pdb_ids_pruned)

df_prot.to_csv(snakemake.output[0], sep='\t', index=False)