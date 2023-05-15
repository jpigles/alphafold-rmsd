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
    df['pdb'] = ''

    for i in range(len(df)):
        # Define UniProt ID and URL
        uniprot_id = df.loc[i, 'uniprot']
        url = 'https://search.rcsb.org/rcsbsearch/v2/query'
    
        pdb_ids = utils.query_rcsb(uniprot_id, url)

        # If received NaN from query, then drop row. Else, prune chains.
        if type(pdb_ids) == float:
            df = df.drop(index=[i])
        else:
            pdb_ids_pruned = utils.prune_extra_chains(pdb_ids)
            df.loc[i, 'pdb'] = pdb_ids_pruned
    
    return df

def download_pdb_files(df, path):
    '''Downloads the PDB files for each protein in the dataframe and saves them in a directory with the Uniprot ID.'''
    for i in range(len(df)):
        uniprot = df.loc[i, 'uniprot']
        uniprot_path = path + uniprot + '/'
        
        # Try to make a new directory with the gene name. If such a directory
        # already exists then continue
        try:
            os.mkdir(uniprot_path)
        except:
            continue
        
        pdb_ids_chains = df.loc[i, 'pdb']
        
        # Remove chains from the PDB IDs
        pdb_ids_no_chains = utils.remove_chains(pdb_ids_chains)

        # A PDB list object that allows to download PDB files
        pdbl = PDBList(verbose=False)

        print('Downloading structures for %s' % uniprot)

        # Retrieve the PDB file from the PDB and save to the directory with the gene name
        pdbl.download_pdb_files(pdb_ids_no_chains, pdir=uniprot_path, file_format='mmCif')

def correct_offset(df, path):

    offsets = []

    for i in range(len(df)):
        # Designate values for retrieval
        uniprot = df.loc[i, 'uniprot']
        pdb_id = df.loc[i, 'pdb']
  
        # Designate file locations. Note that we will be overwriting the CIF files
        cif_path = path + uniprot + '/' + pdb_id + '.cif'

        # Get the offset
        offset = utils.get_offset(cif_path, pdb_id)
        offsets.append(offset)

        # Fix the offset
        fixed_pdb = utils.fix_offset(pdb_id, cif_path, offset)
        print(fixed_pdb)

    # Add column with offset values
    df.insert(len(df.columns), 'auth_offset', offsets)

    return df

def find_domain_completeness(df, path):

    # Convert the domain region strings to ranges or lists of ranges
    df['region_1 search'] = df['region_1'].apply(lambda x: utils.string2range(x))
    df['region_2 search'] = df['region_2'].apply(lambda x: utils.string2range(x))

    for i in range(len(df)):
    
        # Define values for retrieval
        region_1_res = df.loc[i, 'region_1 search']
        region_2_res = df.loc[i, 'region_2 search']
        pdb = df.loc[i, 'pdb']
        uniprot = df.loc[i, 'uniprot']
        path_uniprot = path + uniprot + '/'
        chain = df.loc[i, 'chain']

        # Get structure and dictionary objects
        structure, mmcif_dict = utils.get_structure_dict(pdb, path_uniprot)

        if mmcif_dict['_exptl.method'][0] == 'X-RAY DIFFRACTION':
            resolution = float(mmcif_dict["_refine.ls_d_res_high"][0])
        elif mmcif_dict['_exptl.method'][0] == 'SOLUTION NMR':
            resolution = np.nan

        # Count number of residues in each region
        count_res_reg_1, count_res_reg_2 = utils.count_residues(region_1_res, region_2_res, structure, chain)

        # Calculate the percentage of residues in each region
        percent_reg_1, percent_reg_2 = utils.calculate_domain_completeness(region_1_res, region_2_res, count_res_reg_1, count_res_reg_2)
