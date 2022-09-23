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

    #Designate Uniprot ID for this loop
    uniprot_id = df_uniprot.loc(i, "Uniprot_ID")

    #Designate file path for Alphafold structure
    path_alphafold = 'data/alphafold_files/'

    #make a parser object
    parser = MMCIFParser(QUIET=True)

    #Create the Alphafold structure object
    af_structure = parser.get_structure(uniprot_id, path_alphafold + 'F-' + uniprot_id + '-F1-model_v3.cif')

    #Make an MMCIFDict structure to get more information from it
    mmcif_dict_af = MMCIF2Dict((path_alphafold + 'F-' + uniprot_id + '-F1-model_v3.cif'))

    #Make a cealigner object
    aligner = CEAligner()

    #Get the guide coordinates for the AF 

    #Set the AlphaFold structure as the reference structure
    ref_structure = cealign.set_reference(af_structure)

    #Create the set of PDB files
    pdb_ids = set(df_uniprot.loc(i, "PDB").split(" "))

    #name of the PDB file
    pdb_id = df_prot.loc(i, "PDB ID")
    
    #define the file path of the PDB file
    uniprot = df_prot.loc(i, "Uniprot_ID")
    path_pdb = 'data/structures/%u/' % uniprot

   
    #I don't want to make a new structure object for my Alphafold file for every single PDB file. Any way I can make it once for each unique Uniprot and compare the PDBs against it?


