from biopandas.mmcif import PandasMmcif
import pandas as pd

df = pd.read_csv('../data/proteins_pdb_best.tsv', sep='\t').astype('object')

auth_chains = []

for i in range(len(df)):
    pdb = df.loc[i, 'PDB']
    cif_path = (f'../data/input/RCSB_cif/{pdb}.cif')
    pmcif = PandasMmcif()
    pred = pmcif.read_mmcif(cif_path)