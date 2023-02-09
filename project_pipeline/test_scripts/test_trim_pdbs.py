from biopandas.pdb import PandasPdb
import pandas as pd

pdb_list = pd.read_csv('./data/proteins_pdb_best.tsv', sep='\t').astype('object')

for i in range(len(pdb_list)):
    for i in range(len(df)):
    model = df.loc(i, 'Model')
    if model != 0:
        continue
    uniprot = df.loc(i, 'Uniprot_ID')
    pdb = df.loc(i, 'PDB ID')
    auth_chain = df.loc(i, 'Auth_chain')
    label_chain = df.loc(i, 'Label_chain')
    gt_path = f'./sample_data/input/sample_pdbs/{pdb}.pdb'
    pred_path = f'./data/output/RCSB_af_full/poly_g_6_fasta/{pdb}.fasta/ranked_0.pdb'