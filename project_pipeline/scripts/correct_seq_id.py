from biopandas.pdb import PandasPdb
import pandas as pd
import requests

#initiate Pandas object
ppdb = PandasPdb()

# Load list of pdb files
pdb_list = pd.read_csv('./project_pipeline/data/proteins_pdb_best.tsv', sep = '\t').astype('object')

# Get offset between author and 

for i in range(len(pdb_list)):
    pdb_id = pdb_list.loc(i, 'PDB ID')

    pdb_df = ppdb.read_csv(f'./project_pipeline/data/input/RCSB_best/{pdb_id}.pdb')
