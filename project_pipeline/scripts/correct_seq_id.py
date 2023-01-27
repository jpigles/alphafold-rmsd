from biopandas.pdb import PandasPdb
import pandas as pd
import requests
uniprot = 'Q16644'

#To access our chain of interest, an example splice can be seen below
# print(req_json[uniprot]['mappings'][0]['segments'][0]['chains'])

# Load list of pdb files
pdb_list = pd.read_csv('./project_pipeline/data/proteins_pdb_best.tsv', sep = '\t').astype('object')

# Get offset between author and uniprot seq id
offsets = []
for i in range(len(pdb_list)):

    # Info needed for get request
    ent_id = pdb_list.loc[i, 'Entity_id']
    pdb_id = pdb_list.loc[i, 'PDB ID']
    url = f'https://www.ebi.ac.uk/pdbe/graph-api/pdbe_pages/uniprot_mapping/{pdb_id}/{ent_id}'

    # Send requests
    print(f'Trying {pdb_id}...')
    req = requests.get(url=url)
    print('Status: {status} for PDB {pdb}'.format(status=req.status_code, pdb=pdb_id))
    if req.status_code != 200:
        continue
    
    # Extract uniprot and startIndex ids, get offset (defined as author start - uniprot start)
    req_json = req.json()

    startIndex = req_json[pdb_id]['data'][0]['residues'][0]['startIndex']
    unpStart = req_json[pdb_id]['data'][0]['residues'][0]['unpStartIndex']
    offset = startIndex - unpStart
    offsets = offsets.append(offset)
    print(f'For {pdb_id}, the offset is {offset}')

    #initiate Pandas object
    ppdb = PandasPdb()




    

for i in range(len(pdb_list)):
    pdb_id = pdb_list.loc(i, 'PDB ID')

    pdb_df = ppdb.read_csv(f'./project_pipeline/data/input/RCSB_best/{pdb_id}.pdb')
