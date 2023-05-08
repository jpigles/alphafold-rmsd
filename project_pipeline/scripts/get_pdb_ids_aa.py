# -*- coding: utf-8 -*-
"""
Created on Thu May 14 11:18:27 2020

@author: Jorge Holguin
"""

from Bio.PDB.PDBList import PDBList
import pandas as pd
import numpy as np
import utils
import os

# Define the download path for the CIF files
cif_path = 'data/input/RCSB_cif/'
df_prot = pd.read_csv(snakemake.input[0], sep = '\t')
# df_prot = pd.read_csv('../data/protein_list.tsv', sep = '\t')

# Create a column to store the PDB IDs for each protein
df_prot['PDB'] = ''

# Iterate through the rows of df_prot
for i in range(len(df_prot)):
    
  # Define UniProt ID and URL
  uniprot_id = df_prot.loc[i, 'Uniprot_ID']
  url = 'https://search.rcsb.org/rcsbsearch/v2/query'
    
  pdb_ids = utils.query_rcsb(uniprot_id, url)

  # If received NaN from query, then drop row. Else, prune chains.
  if type(pdb_ids) == float:
    df_prot = df_prot.drop(index=[i])
  else:
    pdb_ids_pruned = utils.prune_extra_chains(pdb_ids)
    df_prot.loc[i, 'PDB'] = pdb_ids_pruned
        
# Save the df_prot as a tsv file
df_prot.to_csv(snakemake.output[0], sep = '\t', index = False)

# Download the pdb files
for i in range(len(df_prot)):
    uniprot = df_prot.loc[i, 'Uniprot_ID']
    uniprot_path = cif_path + uniprot + '/'
    
    # Try to make a new directory with the gene name. If such a directory
    # already exists then continue
    try:
        os.mkdir(uniprot_path)
    except:
        continue
    
    pdb_ids_chains = df_prot.loc[i, 'PDB']
    
    # Remove chains from the PDB IDs
    pdb_ids_no_chains = utils.remove_chains(pdb_ids_chains)

    # A PDB list object that allows to download PDB files
    pdbl = PDBList(verbose=False)

    print('Downloading structures for %s' % uniprot)

    # Retrieve the PDB file from the PDB and save to the directory with the gene name
    pdbl.download_pdb_files(pdb_ids_no_chains, pdir=uniprot_path, file_format='mmCif')