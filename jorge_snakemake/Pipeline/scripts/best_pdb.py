# -*- coding: utf-8 -*-
"""
Created on Wed May 13 16:24:49 2020

@author: Jorge Holguin
"""

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import pandas as pd
import numpy as np
import re
from mutation_enrichment import string2range

# path = 'C:\\Users\\Jorge Holguin\\Documents\\UBC\\4. Fourth Year\\BIOC 448\\Data\\Structures\\Data Files\\'
# path_IAS = '../data/protein_list_pdb.tsv'

# Create an empty dataframe to store the information from the pdb files
df_pdb = pd.DataFrame(columns = ['Gene_name', 'Uniprot_ID', 'Protein_length', 'region_1', 'region_2', 'region_1_len', 
                                 'region_2_len', 'PDB ID', 'PDB Length', 'Resolution',
                                 'Model', 'Chain', 'PDB residues in region_1', 'PDB residues in region_2', 
                                 'Percent residues in region_1', 'Percent residues in region_2'])

# Open the data of proteins with an IAS
df_prot = pd.read_csv(snakemake.input[0], sep = '\t').astype('object')
# df_prot = pd.read_csv('../data/protein_list_pdb.tsv', sep = '\t').astype('object')

# Keep only the rows which have a PDB file
df_prot = df_prot.dropna(subset = ['PDB']).reset_index(drop = True)

# Go through the region_1 and region_2 column and convert the strings into ranges or lists of 
# ranges
df_prot['region_1 search'] = df_prot['region_1'].apply(lambda x: string2range(x))
df_prot['region_2 search'] = df_prot['region_2'].apply(lambda x: string2range(x))    

# Go through the df_prot and each of the PDB files of the proteins with an IAS to 
# determine the number of residues inside the IAS and the Domain that are in the 
# PDB files
for i in range(len(df_prot)):
    
    region_1_res = df_prot.loc[i, 'region_1 search']
    region_2_res = df_prot.loc[i, 'region_2 search']

    # Get the names of all the PDB files corresponding to one protein in a set
    pdb_ids = set(df_prot.loc[i, 'PDB'].split(sep = ','))
    
    # Define the file path for the PDB files of that protein
    gene = df_prot.loc[i, 'Gene_name']
    path_gene = snakemake.input[i+1] + '/'
    # path_gene = '../data/structures/' + gene + '/'
    
    if gene in path_gene:
    
        print('Determining the best structures for %s' % gene)
    
        # Iterate through all the PDB files and determine the number of residues inside
        # the IAS and the Domain that are present in each PDB file
        for pdb in pdb_ids:
            
            # Analyze only the proteins with one unique chain. This line takes the 
            # last characters in the string, retrieves the intergers and checks whether
            # it is equal to 1
            if 1 == int((re.findall('\d+', pdb[4:]))[0]):
        
                # The PDB ids are in the format pdb_id:unique_chains, we only want the id
                pdb = pdb[:4]
                
                # To load a PDB file make a parser object
                parser = MMCIFParser(QUIET=True)
                
                # Then make a structure object
                structure = parser.get_structure(pdb, path_gene + pdb + '.cif')
                
                # Make an MMCIFDict object to grab more information form the .cif files
                mmcif_dict = MMCIF2Dict((path_gene + pdb + '.cif'))
                
                # If the structure was determined through X-Ray crystallography, get the resolution of the structure
                # If the structure was determined through NMR, leave blank
                if mmcif_dict['_exptl.method'][0] == 'X-RAY DIFFRACTION':
                    resolution = float(mmcif_dict["_refine.ls_d_res_high"][0])
                elif mmcif_dict['_exptl.method'][0] == 'SOLUTION NMR':
                    resolution = np.nan
                
                # Iterate through all the models in the structure object (useful only for NMR structures)
                for model in structure:
                    
                    for chain in model:
                    
                        # Get all the residues in the chain A
                        residues = chain.get_residues()
                        
                        # Set all the counters to zero
                        count_res = 0
                        count_res_region_1 = 0
                        count_res_region_2 = 0
                        count_dis_res_region_1 = 0
                        count_dis_res_region_2 = 0
        
                        # Iterate through all the residues in the chain and determine
                        # whether they belong to the IAS or to the Domain. Also determine 
                        # the number of disordered residues in the IAS and in the Domain
                        for residue in residues:
                            count_res = count_res + 1                    
                            
                            # Amino acid residues have an empty space in position zero
                            # of the id
                            if residue.get_id()[0] == ' ':
                                # The sequence position of the amino acid residue is stored
                                # in position 1 of the id
                                if residue.get_id()[1] in region_1_res:
                                    # print(residue.get_id()[1])
                                    count_res_region_1 = count_res_region_1 + 1
                                    
                                    if residue.is_disordered() == 1:
                                        count_dis_res_region_1 = count_dis_res_region_1 + 1
                                    
                                elif residue.get_id()[1] in region_2_res:
                                    count_res_region_2 = count_res_region_2 + 1
                                    
                                    if residue.is_disordered() == 1:
                                        count_dis_res_region_2 = count_dis_res_region_2 + 1
                        
                        # Calculate the percentage of residues in the PDB structure that
                        # occur within the IAS or the Domain
                        percent_in_region_1 = (count_res_region_1/len(region_1_res))*100
                        percent_in_region_2 = (count_res_region_2/len(region_2_res))*100
                        percent_dis_in_region_1 = (count_dis_res_region_1/len(region_1_res))*100
                        percent_dis_in_region_2 = (count_dis_res_region_2/len(region_2_res))*100
                        
                        # Append the results to the df_pdb
                        df_pdb = df_pdb.append({'Gene_name': df_prot.loc[i, 'Gene_name'],
                                                'Uniprot_ID': df_prot.loc[i, 'Uniprot_ID'],
                                                'Protein_length': df_prot.loc[i, 'Protein_length'],
                                                'region_1': df_prot.loc[i, 'region_1'],
                                                'region_2': df_prot.loc[i, 'region_2'],
                                                'region_1_len': len(region_1_res),
                                                'region_2_len': len(region_2_res),
                                                'PDB ID': pdb, 
                                                'PDB Length': count_res, 
                                                'Resolution': resolution, 
                                                'Model': model.get_id(), 
                                                'Chain': chain.get_id(),
                                                'PDB residues in region_1': count_res_region_1,
                                                'PDB residues in region_2': count_res_region_2,
                                                'Percent residues in region_1': percent_in_region_1,
                                                'Percent residues in region_2': percent_in_region_2}, ignore_index=True)
                
# Store the files where more than 80% of the IAS and the Domain exist in the structure
df_pdb_best = df_pdb.loc[(df_pdb['Percent residues in region_1'] > 80.0) & (df_pdb['Percent residues in region_2'] > 80.0)].reset_index(drop = True)

# Save the file with all the pdb files
# df_pdb.to_csv(path + 'pdb_summary.tsv', sep = '\t', index = False)

# Save the file with all the best pdb files
df_pdb_best.to_csv(snakemake.output[0], sep = '\t', index = False)