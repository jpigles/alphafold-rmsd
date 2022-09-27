import pandas as pd

from pdbecif.mmcif_io import CifFileWriter
from biopandas.mmcif import PandasMmcif

#open some .cif file using biopandas mmcif.
pcif = PandasMmcif()
af = pcif.read_mmcif('project_pipeline/data/input/alphafold_files/F-O43739-F1-model_v3.cif')
df = pcif.df['ATOM']
print(df.head(5))
afdict = pd.DataFrame.to_dict(df)