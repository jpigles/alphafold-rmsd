from biopandas.mmcif import PandasMmcif
import pandas as pd

df = pd.read_csv('../data/proteins_pdb_best.tsv', sep='\t').astype('object')

auth_chains = []
