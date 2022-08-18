# -*- coding: utf-8 -*-
"""
Created on Thu May 14 15:00:38 2020

@author: Jorge Holguin
"""

from Bio.PDB.PDBList import PDBList
import pandas as pd
import math
import os

# path = './data/structures/'

df_prot = pd.read_csv(snakemake.input[0], sep = '\t')

for item in range(len(df_prot)):
    gene = df_prot.loc[item, 'Gene_name']
    
    if gene in snakemake.output[item]:
    
        # Try to make a new directory with the gene name. If such a directory
        # already exists then continue
        try:
            os.mkdir(snakemake.output[item])
        except:
            continue
        
        pdb_ids = df_prot.loc[item, 'PDB']
        
        if pd.isna(pdb_ids):
            
            print('No structures found for %s' % gene)
            
        else:
            
            pdb_ids = pdb_ids.split(sep = ',')

            pdb_ids = [i[:-2] for i in pdb_ids]

            # A PDB list object that allows to download PDB files
            pdbl = PDBList(verbose=False)

            print('Downloading structures for %s' % gene)

            # Retrieve the PDB file from the PDB and save to the directory with the gene name
            pdbl.download_pdb_files(pdb_ids, pdir=snakemake.output[item], file_format='mmCif')
