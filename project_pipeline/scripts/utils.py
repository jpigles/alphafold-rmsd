from pdbecif.mmcif_io import CifFileReader, CifFileWriter
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

def remove_chains(pdb_ids_chains):
  #The pdb ids will have their chains attached here (format example: 5ecy.A)
  pdb_ids_chains_list = pdb_ids_chains.split(sep=' ')

  #empty list to store pdb ids without chains
  pdb_ids_no_chains = []

  #Remove the chains from the PDB ids
  for pdb_id in pdb_ids_chains_list:
      pdb_id_only = pdb_id.split(sep='.')[0]
      pdb_ids_no_chains.append(pdb_id_only)
  
  return pdb_ids_no_chains

def expand_on_pdbs(df):
  # Convert PDB column to list-like
  df['pdb'] = df['pdb'].str.split(sep=' ')
  # Explode PDB column
  df = df.explode('pdb').reset_index(drop = True)
  # Split PDB ID and chain into separate columns
  df[['pdb', 'label_chain']] = df['pdb'].str.split(sep='.', expand = True)

  return df

def get_offset(fp, pdb):
  # initiate reader object
  cfr = CifFileReader()
  cif_obj = cfr.read(fp, output='cif_wrapper')
  cif_data = list(cif_obj.values())[0]
  
  # Extract the auth_seq start and db_seq start (from Uniprot) from the cif file
  auth_start = int(cif_data._struct_ref_seq.pdbx_auth_seq_align_beg[0])
  unp_start = int(cif_data._struct_ref_seq.db_align_beg[0])
  offset = auth_start - unp_start

  print(f'Offset for {pdb}: {offset}')
  return offset

def fix_offset(pdb, fp, chain, offset):
    
    # Replace residue numbers in our chain of interest
    if offset == 0:
        return f'No fix needed for {pdb}'
    else:
        # Read in the CIF file
        cfr = CifFileReader()
        cif_obj = cfr.read(fp, output='cif_dictionary', only=['_atom_site'])
        # Convert to a Pandas DataFrame
        df = pd.DataFrame.from_dict(cif_obj[pdb.upper()]['_atom_site'])
        # Replace residue numbers
        for i in range(len(df)):
            if df.loc[i, 'chain_id']==chain:
                res_num = df.loc[i, 'residue_number']
                new_res_num = res_num - offset
                df.loc[i, 'residue_number'] = new_res_num
            else:
                continue
        # Convert back to mmCIF-like dictionary
        cif_dict = df.to_dict(orient='list')
        cif_obj[pdb.upper()]['_atom_site'] = cif_dict

        # Write to file
        cfw = CifFileWriter(fp)
        cfw.write(cif_obj)

        return f'Successfully fixed {pdb}'