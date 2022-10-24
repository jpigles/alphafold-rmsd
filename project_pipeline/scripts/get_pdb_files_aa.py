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
        
        pdb_ids = df_prot.loc[item, 'PDB']
        
        if pd.isna(pdb_ids):
            
            print('No structures found for %s' % uniprot)
            
        else:
            
            #The separator for the PDB IDs I have in my file is a space. I can reconfigure this to make it a comma.
            pdb_ids = pdb_ids.split(sep = ' ')

            # pdb_ids = []
            # for pdb_id_lower in pdb_ids_lower:
            #     pdb_id = pdb_id_lower.upper()
            #     pdb_ids.append(pdb_id)

            # pdb_ids = [i[:-2] for i in pdb_ids] #Does this get rid of the comma?

            # A PDB list object that allows to download PDB files
            pdbl = PDBList(verbose=False)

            print('Downloading structures for %s' % uniprot)

            # Retrieve the PDB file from the PDB and save to the directory with the gene name
            pdbl.download_pdb_files(pdb_ids, pdir=snakemake.output[item], file_format='mmCif')
