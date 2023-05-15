# -*- coding: utf-8 -*-
"""
Created on Thu May 14 11:18:27 2020

@author: Jorge Holguin

Edited by Brooks Perkins-Jechow
"""

from Bio.PDB.PDBList import PDBList
import pandas as pd
import numpy as np
import utils
import main
import os

# Define the download path for the CIF files
cif_path = 'data/input/RCSB_cif/'
df_prot = pd.read_csv(snakemake.input[0], sep = '\t')
# df_prot = pd.read_csv('../data/protein_list.tsv', sep = '\t')

print('Querying RCSB for PDB IDs.')

df_prot = main.get_pdb_ids(df_prot)
print('Successfully retrieved IDs. Proceeding to download structures.')

# Download the pdb files
df_prot = main.download_pdb_files(df_prot, cif_path)

# Make a new dataframe with each PDB ID in a separate row and chains in their own column
df_prot = utils.expand_on_pdbs(df_prot)

# Fix any offsets between the UniProt sequence and the PDB sequence in the CIF files
df_prot = main.correct_offset(df_prot)

# Save the dataframe as a tsv file
df_prot.to_csv(snakemake.output[0], sep = '\t', index = False)

