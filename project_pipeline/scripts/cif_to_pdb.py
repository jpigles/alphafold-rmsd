import pymol
import pandas as pd
from pdb-tools import pdb_fromcif

#We consult the list of uniprot ID's in best_pdb.
#Then we iterate through those in our RCSB_cif and alphafold_files_cif\
#One by one, we turn them into PDB files and insert them into RCSB_pdb and Alphafold_pdb

# prot_df = pd.read_csv(snakemake.input[0], sep = '\t').astype('object')

# for i in len(prot_df.loc('PDB ID')):
#     pdb = prot_df.loc(i, 'PDB ID')
#     uniprot = prot_df.loc(i, 'Uniprot ID')
#     filename = '../data/input/RCSB_cif/{uniprot}/{PDB}.cif'.format(uniprot=uniprot, pdb=pdb)
#     with pymol2.PyMOL() as pymol:
#         pymol.cmd.load(filename)
#         pymol.cmd.save(filename.replace('cif', 'pdb'))

# df_by_uniprot = prot_df.sort_values('Uniprot ID', ascending=False)
# df_unique_uniprots = df.drop_duplicates(subset='Uniprot ID', keep='first')
# print('Sorting by unique Uniprot IDs...')

# for i in len(df_unique_uniprots.loc('Uniprot_ID')):
#     uniprot = df_unique_uniprots.loc(i, 'Uniprot_ID')
#     filename = '../data/input/Alphafold_cif/F-{uniprot}-F1-model_v3.cif'.format(uniprot = uniprot)
#     with pymol2.PyMOL() as pymol:
#         pymol.cmd.load(filename)
#         pymol.cmd.save(filename.replace('cif', 'pdb')) #This will replace both Alphafold_cif with Alphafold_pdb and file.cif with file.pdb
