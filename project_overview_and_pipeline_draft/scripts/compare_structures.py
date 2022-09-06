"""
Created on Sep 6 2022 14:04
@author: Brooks Perkins-Jechow
"""

from Bio.PDB import MMCIFParser, cealign, NeighborSearch, Selection
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import csv
import pandas as pd

#Open the data of proteins with PDB files associated
df_prot = pd.read_csv(snakemake.input[0], sep="\t").astype("object")

#Open the data with the list of PDBs for each Uniprot in the same cell
df_uniprot = pd.read_csv(snakemake.input[1], sep=",").astype("object")

#Loop through each unique Uniprot
for i in range(len(df_uniprot)):

    #Create the list of PDB files

    #Create a parser object

    #Create the Alphafold structure object

    #Set it as the reference structure

    #name of the PDB file
    pdb_id = df_prot.loc(i, "PDB ID")
    
    #define the file path of the PDB file
    uniprot = df_prot.loc(i, "Uniprot_ID")
    path_pdb = 'data/structures/%u/' % uniprot

    #define the file path of the PDB file
    path_alphafold = 'data/alphafold_files/'

    #make a parser object
    parser = MMCIFParser(QUIET=True)

    #I don't want to make a new structure object for my Alphafold file for every single PDB file. Any way I can make it once for each unique Uniprot and compare the PDBs against it?


