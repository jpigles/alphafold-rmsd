# -*- coding: utf-8 -*-
"""
Created on Thu May 14 11:18:27 2020

@author: Jorge Holguin
"""
 
import requests
import pandas as pd
import numpy as np

df_prot = pd.read_csv(snakemake.input[0], sep = '\t')
# df_prot = pd.read_csv('../data/protein_list.tsv', sep = '\t')

# Create a column to store the PDB IDs for each protein
df_prot['PDB'] = ''

# Iterate through the rows of df_prot
for i in range(len(df_prot)):
    
    # Create a variable with the Uniprot_ID
    uniprot_id = df_prot.loc[i, 'Uniprot_ID']

    url = 'http://www.rcsb.org/pdb/rest/search'
    
    # Query text to request all the PDB IDs associated to a Uniprot_ID
    # Taken from https://www.rcsb.org/pdb/software/rest.do
    query_text = """
<?xml version="1.0" encoding="UTF-8"?>

<orgPdbQuery>

<queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>

<description>Simple query for a list of Uniprot Accession IDs</description>

<accessionIdList>%s</accessionIdList>

</orgPdbQuery>

""" % uniprot_id
    
    print("Querying RCSB PDB REST API for Uniprot_ID: %s" % uniprot_id)
    
    header = {'Content-Type': 'application/x-www-form-urlencoded'}
    
    response = requests.post(url, data=query_text, headers=header)
    
    # Store the results in df_prot
    if response.status_code == 200:
        pdb_id = response.text.strip().replace('\n', ',')
        df_prot.loc[i, 'PDB'] = pdb_id
        # print("Found %d PDB entries matching query." % len(response.text))
        # print("Matches: \n%s" % response.text)
        
    else:
        pdb_id = np.nan
        df_prot.loc[i, 'PDB'] = pdb_id
        print("Failed to retrieve results for Uniprot_ID: %s" % uniprot_id)
        
# Save the df_prot as a tsv file
df_prot.to_csv(snakemake.output[0], sep = '\t', index = False)
        
