'''
Created on Friday Oct 14th 2022

Author: Brooks Perkins-Jechow

'''

import pandas as pd

def prune_extra_chains(pdb_ids_str):

    #Turn the string into a list
    pdb_ids_w_chain = pdb_ids_str.strip().split(sep=' ')

    #Empty dictionary to fill.
    pdb_ids_dict = {}

    for pdb_id in pdb_ids_w_chain:

        #PDB ID (lowercase)
        pdb = pdb_id[:4].lower()

        #Chain label
        chain = pdb_id[5]

        #Add the PDB ID as a key and the chain label as a value.
        if pdb not in pdb_ids_dict.keys():
            pdb_ids_dict[pdb] = [chain]
        else:
            pdb_ids_dict[pdb].append(chain)

    #Extract each PDB from the pdb_ids dict
    for pdb_id in pdb_ids_dict.copy():
            
        #Determine whether there are one or more chains. If there are more than one chain, pop it out.
        if len(pdb_ids_dict[pdb_id]) != 1:

            #Remove the PDB ID from the pdb_ids_dict
            remove_pdb = pdb_ids_dict.pop(pdb_id)

    # Now we convert them back to strings and add it all together.
    # Make a list of the values
    values_list = list(pdb_ids_dict.values())

    #Make a list of the keys
    key_list = list(pdb_ids_dict.keys())

    #Empty string to fill with my IDs.
    unique_pdb_ids = ''

    for n in range(len(pdb_ids_dict)):

        #Get the value
        chain = values_list[n][0]

        #Get the key
        key = key_list[n].lower()

        # Put them together
        pdb_id_chain_str = key + '.' + chain

         #Append to my unique_pdb_ids string
        unique_pdb_ids = unique_pdb_ids + ' ' + pdb_id_chain_str

    #Make the value of PDB at the index i equal to our new string.
    return unique_pdb_ids


df_prot = pd.read_csv(snakemake.input[0], sep='\t').astype('object')

#Curate the list of PDB IDs
for i in range(len(df_prot)):
    if type(df_prot.loc[i, 'PDB']) == float:
        df_prot = df_prot.drop(index=[i])
    else:
        pdb_ids_unpruned = df_prot.loc[i, 'PDB']
        pdb_ids_pruned = prune_extra_chains(pdb_ids_unpruned)
        df_prot.loc[i, 'PDB'] = pdb_ids_pruned

df_prot.to_csv(snakemake.output[0], sep='\t', index=False)