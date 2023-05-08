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


def prune_extra_chains(pdb_ids_str):

    #Turn the string into a list
    pdb_ids_w_chain = pdb_ids_str.strip().split(sep=' ')

    #Empty dictionary to fill.
    pdb_ids_dict = {}

    for pdb_id in pdb_ids_w_chain:

        #PDB ID (lowercase)
        pdb = pdb_id[:4].lower()

        #Chain label
        chain = pdb_id[5]

        #Add the PDB ID as a key and the chain label as a value.
        if pdb not in pdb_ids_dict.keys():
            pdb_ids_dict[pdb] = [chain]
        else:
            pdb_ids_dict[pdb].append(chain)

    #Extract each PDB from the pdb_ids dict
    for pdb_id in pdb_ids_dict.copy():
            
        #Determine whether there are one or more chains.
        if len(pdb_ids_dict[pdb_id]) != 1:

            #If more than one chain, select the first chain as our representative.
            pdb_ids_dict[pdb_id] = pdb_ids_dict[pdb_id][0]

    # Now we convert them back to strings and add it all together.
    # Make a list of the values
    values_list = list(pdb_ids_dict.values())

    #Make a list of the keys
    key_list = list(pdb_ids_dict.keys())

    #Empty string to fill with my IDs.
    unique_pdb_ids = ''

    for n in range(len(pdb_ids_dict)):

        #Get the value
        chain = values_list[n][0]

        #Get the key
        key = key_list[n].lower()

        # Put them together
        pdb_id_chain_str = key + '.' + chain

         #Append to my unique_pdb_ids string
        unique_pdb_ids = unique_pdb_ids + ' ' + pdb_id_chain_str

    #Make the value of PDB at the index i equal to our new string.
    return unique_pdb_ids