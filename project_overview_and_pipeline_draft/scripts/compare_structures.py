"""
Created on Sep 6 2022 14:04
@author: Brooks Perkins-Jechow
"""

from Bio.PDB import MMCIFParser, cealign, NeighborSearch, Selection
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import csv
import pandas as pd

#Open the data of proteins with PDB files associated
df_uniprot = pd.read_csv(snakemake.input[0], sep=",").astype("object")
