import requests
import pandas as pd
import numpy as np

def query_rcsb(uniprot_id, url):
    
    # Query test for pdb files associated with given UniProt accession number.
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

    # In format 4CJ0 1P3I ...
    pdb_str = ''

    if response.status_code == 200:
        response_dic = response.json()
        for n in range(len(response_dic['result_set'])):
            pdb_str = pdb_str + response_dic['result_set'][n]['identifier'] + ' '
        
    else:
        pdb_str = np.nan
        print("Failed to retrieve results for Uniprot_ID: %s" % uniprot_id)

    return pdb_str