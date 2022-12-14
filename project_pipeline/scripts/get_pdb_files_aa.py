# -*- coding: utf-8 -*-
"""
Created on Thu May 14 15:00:38 2020

@author: Jorge Holguin

Copy created on Tue Aug 23 2022

@author: Brooks Perkins-Jechow
"""

from Bio.PDB.PDBList import PDBList
import pandas as pd
import math
import os

# path = './data/structures/'

df_prot = pd.read_csv(snakemake.input[0], sep = '\t').astype('object')

for item in range(len(df_prot)):
    uniprot = df_prot.loc[item, 'Uniprot_ID']
    
    if uniprot in snakemake.output[item]:
    
        # Try to make a new directory with the gene name. If such a directory
        # already exists then continue
        try:
            os.mkdir(snakemake.output[item])
        except:
            continue
        
        pdb_ids_chains = df_prot.loc[item, 'PDB']
        
        if pd.isna(pdb_ids_chains):
            
            print('No structures found for %s' % uniprot)
            
        else:
            
            #The pdb ids will have their chains attached here (format example: 5ecy.A)
            pdb_ids_chains_list = pdb_ids_chains.split(sep = ' ')

            #empty list to store pdb ids without chains
            pdb_ids_no_chains = []

            #Remove the chains from the PDB ids
            for pdb_id in pdb_ids_chains_list:
                pdb_id_only = pdb_id[:4]
                pdb_ids_no_chains.append(pdb_id_only)

            # pdb_ids = [i[:-2] for i in pdb_ids] #Does this get rid of the comma?

            # A PDB list object that allows to download PDB files
            pdbl = PDBList(verbose=False)

            print('Downloading structures for %s' % uniprot)

            # Retrieve the PDB file from the PDB and save to the directory with the gene name
            pdbl.download_pdb_files(pdb_ids_no_chains pdir=snakemake.output[item], file_format='mmCif')
