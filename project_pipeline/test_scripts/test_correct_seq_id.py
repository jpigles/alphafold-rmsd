from biopandas.pdb import PandasPdb
import pandas as pd
import requests

#To access our chain of interest, an example splice can be seen below
# print(req_json[uniprot]['mappings'][0]['segments'][0]['chains'])

# Load list of pdb files
pdb_list = pd.read_csv('./sample_data/sample_proteins_pdb_best.tsv', sep = '\t').astype('object')

def pdbe_req(ent_id, pdb_id):
    # Send request for pdb id
    url = f'https://www.ebi.ac.uk/pdbe/graph-api/pdbe_pages/uniprot_mapping/{pdb_id}/{ent_id}'
    print(f'Trying {pdb_id}...')
    req = requests.get(url=url)
    print('Status: {status} for PDB {pdb}'.format(status=req.status_code, pdb=pdb_id))
    return req.status_code, req.json()

def get_offset(json, pdb, uniprot):
    # The Uniprot sequence is always taken from the Uniprot-defined canonical sequence (first isoform)
    unpStart = json[pdb]['Uniprot'][uniprot]['mappings'][0]['unp_start']
    unpEnd = json[pdb]['Uniprot'][uniprot]['mappings'][0]['unp_end']
    # There's no way to know which author_provided number will be available, if any, so
    # we try both isoforms
    try:
        auth_start = json[pdb]['Uniprot'][uniprot]['mappings'][0]['start']['author_residue_number']
        offset = auth_start - unpStart
        print('Auth_start was successful')
    except TypeError:
        try:
            auth_end = json[pdb]['Uniprot'][uniprot]['mappings'][0]['end']['author_residue_number']
            offset = auth_end - unpEnd
            print('Auth_end was successful')
        except TypeError:
            try:
                auth_start_2 = json[pdb]['Uniprot'][uniprot + '-2']['mappings'][0]['start']['author_residue_number']
                offset = auth_start - unpStart
                print('Auth_start_2 was successful')
            except TypeError:
                try:
                    auth_end_2 = json[pdb]['Uniprot'][uniprot + '-2']['mappings'][0]['end']['author_residue_number']
                    offset = auth_end - unpEnd
                    print('Auth_end_2 was successful')
                except TypeError:
                    print('No available author numbers')
                    offset = 'Null'
                    null_pdbs = null_pdbs.append(pdb)

    print(f'For {pdb_id}, the offset is {offset}')
    return offset

def fix_seq_id(pdb, in_fp, out_fp, chain, offset):
    # initiate Pandas object
    ppdb = PandasPdb()
    _ = ppdb.read_pdb(in_fp)
    pred = ppdb.df['ATOM']

    # Replace residue numbers in our chain of interest
    for i in range(len(pred)):
        if pred.loc[i, 'chain_id']==chain:
            res_num = pred.loc[i, 'residue_number']
            new_res_num = res_num - offset
            pred.loc[i, 'residue_number'] = new_res_num
        else:
            continue
    
    pred.to_pdb(path=out_fp, records=None, gz=False, append_newline=True)

    return f'Successfully fixed {pdb}'

# Get offset between author and uniprot seq id
offsets = []

# List of pdbs with no available author residue numbers via PDBe
null_pdbs = []
for i in range(len(pdb_list)):

    # Info needed for get request
    ent_id = pdb_list.loc[i, 'Entity_id']
    pdb_id = pdb_list.loc[i, 'PDB ID']
    uniprot_id = pdb_list.loc[i, 'Uniprot_ID']
    
    # Info needed for pdb object
    path = f'./data/input/RCSB/pdbs/{pdb_id}.pdb'
    out_path = f'.sample_data/input/sample_pdbs/{pdb_id}.pdb'
    auth_chain = pdb_list.loc[i, 'Auth_chain']

    # send request
    req_status, req_json = pdbe_req(ent_id, pdb_id)
    if req_status != 200:
        continue

    # Get offset
    offset = get_offset(req_json, pdb_id, uniprot_id)
    offsets.append(offset)

    # fix pdb sequence ids
    fixed_pdb = fix_seq_id(pdb_id, path, out_path, auth_chain, offset)
    print(fixed_pdb)

# Add offsets to file
pdb_list.insert(18, 'Auth_offset', offsets)
pdb_list.to_csv('./sample_data/sample_proteins_pdb_best.tsv', sep = '\t', index=False)