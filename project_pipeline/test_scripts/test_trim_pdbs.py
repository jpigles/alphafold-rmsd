from biopandas.pdb import PandasPdb
import pandas as pd

pdb_list = pd.read_csv('./sample_data/sample_proteins_pdb_best.tsv', sep='\t').astype('object')

for i in range(len(pdb_list)):
    # skip extra rows for NMR files
    model = pdb_list.loc(i, 'Model')
    if model != 0:
        continue
    
    # Define parameters for selecting files
    uniprot = pdb_list.loc(i, 'Uniprot_ID')
    pdb = pdb_list.loc(i, 'PDB ID')
    auth_chain = pdb_list.loc(i, 'Auth_chain')
    label_chain = pdb_list.loc(i, 'Label_chain')
    gt_path = f'./sample_data/input/sample_pdbs/{pdb}.pdb'
    pred_path = f'./data/output/RCSB_af_full/poly_g_6_fasta/{pdb}.fasta/ranked_0.pdb'

    # gt model
    ppdb = PandasPdb()
    gt_file = ppdb.read_pdb(gt_path)
    gt = gt_file['ATOM']

    pred_file = ppdb.read_pdb(pred_path)
    pred = pred_file['ATOM']

    