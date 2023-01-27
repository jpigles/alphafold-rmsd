from biopandas.pdb import PandasPdb
import pandas as pd
import requests
uniprot = 'Q16644'

#To access our chain of interest, an example splice can be seen below
# print(req_json[uniprot]['mappings'][0]['segments'][0]['chains'])

# Load list of pdb files
pdb_list = pd.read_csv('./project_pipeline/data/proteins_pdb_best.tsv', sep = '\t').astype('object')

def pdbe_req(ent_id, pdb_id):
    # Send request for pdb id
    url = f'https://www.ebi.ac.uk/pdbe/graph-api/pdbe_pages/uniprot_mapping/{pdb_id}/{ent_id}'
    print(f'Trying {pdb_id}...')
    req = requests.get(url=url)
    print('Status: {status} for PDB {pdb}'.format(status=req.status_code, pdb=pdb_id))
    return req.status_code, req.json()

def get_offset(json):
    startIndex = json[pdb_id]['data'][0]['residues'][0]['startIndex']
    unpStart = json[pdb_id]['data'][0]['residues'][0]['unpStartIndex']
    offset = startIndex - unpStart
    print(f'For {pdb_id}, the offset is {offset}')
    return offset

def fix_seq_id(fp, chain, offset):
    #initiate Pandas object
    ppdb = PandasPdb()
    _ = ppdb.read_pdb(fp)
    pred = ppdb.df['ATOM']

    # Select only our desired chain
    

    for i in range(len(pred)):
        if pred.loc[i, 'chain_id']==chain:
            res_num = pred.loc[i, 'residue_number']
            new_res_num = res_num - offset
            pred.loc[i, 'residue_number'] = new_res_num
        else:
            continue

# Get offset between author and uniprot seq id
offsets = []
for i in range(len(pdb_list)):

    # Info needed for get request
    ent_id = pdb_list.loc[i, 'Entity_id']
    pdb_id = pdb_list.loc[i, 'PDB ID']
    
    # Info needed for pdb object
    path = f'./data/input/RCSB/pdbs/{pdb_id}.pdb'
    auth_chain = pdb_list.loc[i, 'Auth_chain']

    # send request
    req_status, req_json = pdbe_req(ent_id, pdb_id)
    if req_status != 200:
        continue

    # Get offset
    offset = get_offset(req_json)
    offsets.append(offset)

    # fix pdb sequence ids
    fixed_pdb = fix_seq_id(path, auth_chain, offset)

for i in range(len(pdb_list)):
    pdb_id = pdb_list.loc(i, 'PDB ID')

    pdb_df = ppdb.read_csv(f'./project_pipeline/data/input/RCSB_best/{pdb_id}.pdb')
