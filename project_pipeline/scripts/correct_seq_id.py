from biopandas.pdb import PandasPdb
import pandas as pd
import requests
uniprot = 'Q16644'

url = f'https://www.ebi.ac.uk/pdbe/graph-api/uniprot/{uniprot}'

req = requests.get(url=url)

print(req.status_code)

req_json = req.json()

#To access our chain of interest, an example splice can be seen below

print(req_json[uniprot]['mappings'][0]['segments'][0]['chains'])
#initiate Pandas object
ppdb = PandasPdb()

# Load list of pdb files
pdb_list = pd.read_csv('./project_pipeline/data/proteins_pdb_best.tsv', sep = '\t').astype('object')

# Get offset between author and uniprot seq id

for i in range(len(pdb_list)):


for i in range(len(pdb_list)):
    pdb_id = pdb_list.loc(i, 'PDB ID')

    pdb_df = ppdb.read_csv(f'./project_pipeline/data/input/RCSB_best/{pdb_id}.pdb')
