from Bio.PDB.PDBList import PDBList
import pandas as pd
import numpy as np
import utils
import os


###########
# Rule pdb_ids functions

def get_pdb_ids(df):
    '''Retrieves PDB IDs for each protein in the dataframe in the form of ID.chain (e.g. 1A2K.A)'''

    # Create a column to store the PDB IDs for each protein
    df['PDB'] = ''

    for i in range(len(df)):
        # Define UniProt ID and URL
        uniprot_id = df.loc[i, 'Uniprot_ID']
        url = 'https://search.rcsb.org/rcsbsearch/v2/query'
    
        pdb_ids = utils.query_rcsb(uniprot_id, url)

        # If received NaN from query, then drop row. Else, prune chains.
        if type(pdb_ids) == float:
            df = df.drop(index=[i])
        else:
            pdb_ids_pruned = utils.prune_extra_chains(pdb_ids)
            df.loc[i, 'PDB'] = pdb_ids_pruned
    
    return df

def download_pdb_files(df, path):
    for i in range(len(df)):
        uniprot = df.loc[i, 'Uniprot_ID']
        uniprot_path = path + uniprot + '/'
        
        # Try to make a new directory with the gene name. If such a directory
        # already exists then continue
        try:
            os.mkdir(uniprot_path)
        except:
            continue
        
        pdb_ids_chains = df.loc[i, 'PDB']
        
        # Remove chains from the PDB IDs
        pdb_ids_no_chains = utils.remove_chains(pdb_ids_chains)

        # A PDB list object that allows to download PDB files
        pdbl = PDBList(verbose=False)

        print('Downloading structures for %s' % uniprot)

        # Retrieve the PDB file from the PDB and save to the directory with the gene name
        pdbl.download_pdb_files(pdb_ids_no_chains, pdir=uniprot_path, file_format='mmCif')

def correct_offset(df, path):
    