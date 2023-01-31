from biopandas.pdb import PandasPdb
import pandas as pd
import requests
import csv

#To access our chain of interest, an example splice can be seen below
# print(req_json[uniprot]['mappings'][0]['segments'][0]['chains'])

# Load list of pdb files
pdb_list = pd.read_csv('./sample_data/sample_proteins_pdb_best.tsv', sep = '\t').astype('object')


def fix_seq_id(pdb, in_fp, out_fp, chain, offset):
    # initiate Pandas object
    ppdb = PandasPdb()
    _ = ppdb.read_pdb(in_fp)
    pred = ppdb.df['ATOM']

    # Replace residue numbers in our chain of interest
    if offset == 'Null':
        return 'No offset to fix with'
    else:
        for i in range(len(pred)):
            if pred.loc[i, 'chain_id']==chain:
                res_num = pred.loc[i, 'residue_number']
                new_res_num = res_num - offset
                pred.loc[i, 'residue_number'] = new_res_num
            else:
                continue
    
    ppdb.to_pdb(path=out_fp, records=None, gz=False, append_newline=True)

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
    req_status, req_json = pdbe_req(pdb_id)
    if req_status != 200:
        continue

    # Get offset
    offset, nulls = get_offset(req_json, pdb_id, uniprot_id)
    offsets.append(offset)
    null_pdbs = null_pdbs + nulls

    # fix pdb sequence ids
    fixed_pdb = fix_seq_id(pdb_id, path, out_path, auth_chain, offset)
    print(fixed_pdb)

# Add offsets to file
pdb_list.insert(18, 'Auth_offset', offsets)
pdb_list.to_csv('./sample_data/sample_proteins_offset.tsv', sep = '\t', index=False)

with open('null_offset_pdbs.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    for pdb in null_pdbs:
        writer.writerow(pdb)