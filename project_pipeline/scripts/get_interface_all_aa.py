# -*- coding: utf-8 -*-
"""
Created on Thu May 21 10:31:23 2020

@author: Jorge Holguin

Copy created on Wed Aug 24 12:30 2022

@author: Brooks Perkins-Jechow
"""

from Bio.PDB import MMCIFParser, NeighborSearch, Selection
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import pandas as pd
import numpy as np
from mutation_enrichment import string2range

# Open the data of proteins with an region_1
df_prot = pd.read_csv(snakemake.input[0], sep = '\t').astype('object')

# Keep only the rows which have a PDB file
df_prot = df_prot.dropna(subset = ['PDB ID']).reset_index(drop = True)

# Go through the region_1 and region_2 column and convert the strings into ranges or lists of 
# ranges
df_prot['region_1 search'] = df_prot['region_1'].apply(lambda x: string2range(x))
df_prot['region_2 search'] = df_prot['region_2'].apply(lambda x: string2range(x))    

# Create a column to store the interacting residue pairs and another to store the
# residues at the interface
df_prot['PDB Mutations'] = ''
df_prot['Interacting residue pairs'] = ''
df_prot['Interface Residues'] = ''
df_prot['Number Interface Residues'] = ''

# Go through the df_prot and each of the PDB files of the proteins with an region_1 to 
# determine the residues at the interface between the region_1 and the region_2
for i in range(len(df_prot)):
    
    region_1_res = df_prot.loc[i, 'region_1 search']
    region_2_res = df_prot.loc[i, 'region_2 search']
            
    # Get the name of the PDB file 
    pdb_id = df_prot.loc[i, 'PDB ID']
    
    # Define the file path for the PDB file of that protein
    uniprot = df_prot.loc[i, 'Uniprot_ID']
    path_uniprot = 'data/input/RCSB_cif/%s/' % uniprot
    
    # To load a PDB file make a parser object
    parser = MMCIFParser(QUIET=True)
    
    # Then make a structure object
    structure = parser.get_structure(pdb_id, path_uniprot + pdb_id + '.cif')
    
    # Make an MMCIFDict object to grab more information from the .cif files
    mmcif_dict = MMCIF2Dict((path_uniprot + pdb_id + '.cif'))
    
    # Retrieve the mutation information from the .cif file
    df_prot.loc[i, 'PDB Mutations'] = mmcif_dict['_entity.pdbx_mutation'][0]
    
    # Iterate through all the models in the structure
    for model in structure:
        
        # Analyze only the model that corresponds to the current row in df_prot
        if model.get_id() == df_prot.loc[i, 'Model']:
        
            for chain in model:
                
                # Analyze only the chain that corresponds to the current row in df_prot
                if chain.get_id() == df_prot.loc[i, 'Chain']:
                
                    # Get all the atoms in the chain
                    atom_list = Selection.unfold_entities(chain, "A")
                    
                    # Make an empty list to store the atoms in the region_1 and in region_2
                    atoms_ns = []
                    
                    # Iterate through all the atoms in the chain and append the ones that occur inside
                    # the region_1 or the region_2 to the list
                    for atom in atom_list:
                        
                        # Get the parent residue for the atom
                        res = atom.get_parent()
                        
                        # Make sure the residue is an amino acid
                        if res.get_id()[0] == ' ':
                            
                            # Check whether the residue lies inside the region_1 or the region_2 and append
                            # the atoms of these residues into a list
                            if res.get_id()[1] in region_1_res:
                                atoms_ns.append(atom)
                                
                            elif res.get_id()[1] in region_2_res:
                                atoms_ns.append(atom)
                    
                    # Make an NeighborSearch object with all the atoms inside the region_1 and the region_2        
                    ns = NeighborSearch(atoms_ns)
                    
                    # Search for all the interacting residues in the region_1 and in the region_2
                    # with atoms that are within a 6.5 A radius 
                    ns_all = ns.search_all(6.5, 'R')
                    
                    # Make a set to store the residues at the interface
                    interface_res = set()
                    
                    # Save the interacting residue pairs as a list of tuples
                    interacting_pairs = []
                    
                    # Iterate thorugh all the interacting residue pairs and append those that have a residue
                    # in the region_1 and another in the region_2 to a list. Save the residue positions in a set
                    for pairs in ns_all:
                        
                        res_0 = pairs[0].get_id()[1]
                        res_1 = pairs[1].get_id()[1]
                        
                        if res_0 in region_1_res and res_1 in region_2_res:
                            interface_res.add(res_0)
                            interface_res.add(res_1)
                            interacting_pairs.append((res_0, res_1))
                            
                        elif res_1 in region_1_res and res_0 in region_2_res:
                            interface_res.add(res_0)
                            interface_res.add(res_1)
                            interacting_pairs.append((res_0, res_1))
                            
                    # Save the results in the appropriate columns of df_prot
                    if len(interface_res) > 0 and len(interacting_pairs) > 0:    
                        df_prot.loc[i, 'Interacting residue pairs'] = str(interacting_pairs)
                        df_prot.loc[i, 'Interface Residues'] = str(interface_res)
                        df_prot.loc[i, 'Number Interface Residues'] = len(interface_res)
                        
                    else: 
                        df_prot.loc[i, 'Interacting residue pairs'] = np.nan
                        df_prot.loc[i, 'Interface Residues'] = np.nan
                        df_prot.loc[i, 'Number Interface Residues'] = np.nan
                        
# Go through all the proteins in df_prot and determine the pdb file with the greatest
# number of interface residues between the region_1 and the region_2

# Convert the int values to numpy.int64 (this is required for the .idxmax method to be applied)
df_prot.loc[:, 'Number Interface Residues'] = pd.to_numeric(df_prot['Number Interface Residues'])

# Make a new column to flag the rows to keep
df_prot['Keep'] = ''

# Get all the Uniprot_IDs
proteins = set(df_prot['Uniprot_ID'])

# Iterate through all the proteins
for protein in proteins:
    
    print('Determining the interface residues for', protein)
    
    # Grab all the instances of each protein in df_prot
    df_temp = df_prot.loc[df_prot['Uniprot_ID'] == protein]
    
    # Grab the intances of each protein in df_prot with no mutations
    df_temp_no_mut = df_temp.loc[df_temp['PDB Mutations'] == '?']
    
    # Grab all the intances of each protein with mutations
    df_temp_mut = df_temp.loc[df_temp['PDB Mutations'] != '?']
    
    # From the PDB entries with no mutations, select the one with the largest 
    # number of residues at the interface
    if len(df_temp_no_mut) > 0 :
        # Get the index of the max value in the Number Interface Residues column
        max_index_no_mut = df_temp_no_mut['Number Interface Residues'].idxmax(skipna = True)
        
        # Set the row with the max index to True in the Keep column of df_prot
        df_prot.loc[max_index_no_mut, 'Keep'] = True
    
    # If all the PDB entries have mutations, select the one with the largest
    # number of residues at the interface
    elif len(df_temp_mut) > 0:
        # Get the index of the max value in the Number Interface Residues column
        max_index_mut = df_temp_mut['Number Interface Residues'].idxmax(skipna = True)
        
        # Set the row with the max index to True in the Keep column of df_prot
        df_prot.loc[max_index_mut, 'Keep'] = True

# Make a new df with the proteins with the highest number of interface residues
df_prot_keep = df_prot.loc[df_prot['Keep'] == True]
    
# Convert the set or list of residues at the interface into a string
def to_string(x):
    x = str(x)
    
    if x == 'nan':
        return np.nan
    elif '{' in x:
        x = x.replace('{', '').replace('}', '').replace(' ', '')
        return x
    elif '[' in x:
        x = x.replace('[', '').replace(']', '').replace(' ', '')
        return x

df_prot_result = df_prot.copy()
df_prot_keep_result = df_prot_keep.copy()

df_prot_result.loc[:,'Interface Residues'] = df_prot_result['Interface Residues'].apply(to_string)    

df_prot_keep_result.loc[:, 'Interface Residues'] = df_prot_keep_result['Interface Residues'].apply(to_string)

df_prot_keep_result.loc[:, 'Interacting residue pairs'] = df_prot_keep_result['Interacting residue pairs'].apply(to_string)

df_prot_keep_result = df_prot_keep_result.dropna(subset = ['Uniprot_ID']).reset_index(drop = True)
df_prot_keep_result = df_prot_keep_result.drop(['region_1 search', 'region_2 search', 'Keep'], axis = 'columns')

# Save the file with all the best pdb files
df_prot_keep_result.to_csv(snakemake.output[0], sep = '\t', index = False)
                


































