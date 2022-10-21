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

    url = 'https://search.rcsb.org/rcsbsearch/v2/query'
    
    # Query text to request all the PDB IDs associated to a Uniprot_ID
    # Taken from https://www.rcsb.org/pdb/software/rest.do
    query_text = {
  "query": {
    "type": "group",
    "logical_operator": "and",
    "nodes": [
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "exact_match",
          "value": uniprot_id,
          "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession"
        }
      },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "exact_match",
          "value": "UniProt",
          "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name"
        }
      }
    ]
  },
  "request_options": {
    "return_all_hits": True
  },
  "return_type": "polymer_instance"
    }
    
    print("Querying RCSB PDB REST API for Uniprot_ID: %s" % uniprot_id)
    
    header = {'Content-Type': 'application/x-www-form-urlencoded'}
    
    response = requests.post(url, json=query_text)
    
    # Store the results in df_prot
    if response.status_code == 200:
        # pdb_id = response.text.strip().replace('\n', ',')
        # df_prot.loc[i, 'PDB'] = pdb_id
        response_dic = response.json()
        pdb_str = ''
        for n in range(len(response_dic['result_set'])):
            pdb_str = pdb_str + response_dic['result_set'][n]['identifier'] + ' '
        df_prot.loc[i, 'PDB'] = pdb_str
        # print("Found %d PDB entries matching query." % len(response.text))
        # print("Matches: \n%s" % response.text)
        
    else:
        pdb_id = np.nan
        df_prot.loc[i, 'PDB'] = pdb_id
        print("Failed to retrieve results for Uniprot_ID: %s" % uniprot_id)
        
# Save the df_prot as a tsv file
df_prot.to_csv(snakemake.output[0], sep = '\t', index = False)
        
