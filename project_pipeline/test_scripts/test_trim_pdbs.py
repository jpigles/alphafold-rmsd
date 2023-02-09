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
    gt_all_chains = gt_file['ATOM']
    gt = gt_all_chains[gt_all_chains['chain_id'] == auth_chain]

    pred_file = ppdb.read_pdb(pred_path)
    pred = pred_file['ATOM']

    # Create list of rows that are present in both gt and pred based on row indices of pred
    present_atoms = []
    for i in range(len(gt)):
        # Define minimum  parameters to select unique rows
        gt_atom_name = gt.loc[i, 'atom_name']
        gt_chain_id = gt.loc[i, 'chain_id']
        gt_residue_number = gt.loc[i, 'residue_number']

        # Look for matching row in pred
        pred_row = pred.loc[(pred['atom_name'] == gt_atom_name) & (pred['residue_number'] == gt_residue_number) & (pred['chain_id'] == gt_chain_id)]
        if pred_row.empty != True:
            present_atoms.append(pred_row.index)

    # Select all rows in pred not present in gt
    total_atoms = list(pred.index)
    na_atoms = list(set(total_atoms).difference(present_atoms))

    # Create new pred data frame exclusively with atoms present in gt
    pred_trim = pred.drop(index=na_atoms)

    assert len(pred_trim) == gt

        # # pred_chain = pred.loc[pred['chain_id'] == gt_chain_id]
# # pred_residue = pred_chain.loc[pred_chain['residue_number'] == gt_residue_number]
# # pred_atom = pred_residue.loc[pred_residue['atom_name'] == gt_atom_name]
# pred_trim = pred.loc[(pred['chain_id'] == gt_chain_id) & (pred['residue_number'] == gt_residue_number) & (pred['atom_name'] == gt_atom_name)]