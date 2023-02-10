from biopandas.pdb import PandasPdb
import pandas as pd
import numpy as np

pdb_list = pd.read_csv('./sample_data/sample_proteins_pdb_best.tsv', sep='\t').astype('object')


for i in range(len(pdb_list)):
    # skip extra rows for NMR files
    model = pdb_list.loc[i, 'Model']
    if model != 0:
        continue
    
    # Define parameters for selecting files
    uniprot = pdb_list.loc[i, 'Uniprot_ID']
    pdb = pdb_list.loc[i, 'PDB ID']
    auth_chain = pdb_list.loc[i, 'Auth_chain']
    label_chain = pdb_list.loc[i, 'Label_chain']
    gt_path = f'./sample_data/input/sample_pdbs/{pdb}.pdb'
    pred_path = f'./data/output/RCSB_af_full/poly_g_6_fasta/{pdb}.fasta/ranked_0.pdb'
    gt_out_path = f'./sample_data/input/sample_pdbs_trim/{pdb}.pdb'
    pred_out_path = f'./sample_data/output/ds1_af_full/poly_g_20_fasta/{pdb}.fasta/ranked_0.pdb'

    print(f'Trying {pdb}...')

    # gt model
    ppdb = PandasPdb()
    gt_file = ppdb.read_pdb(gt_path)
    gt_all_chains = gt_file.df['ATOM']
    gt = gt_all_chains[gt_all_chains['chain_id'] == auth_chain]

    pred_file = ppdb.read_pdb(pred_path)
    pred = pred_file.df['ATOM']

    print(len(gt), len(pred))

    # Create list of rows that are present in both gt and pred based on row indices of pred
    present_atoms = []
    for atom in range(len(gt)):
        # Define minimum  parameters to select unique rows
        gt_atom_name = gt.loc[atom, 'atom_name']
        gt_residue_name = gt.loc[atom, 'residue_name']
        gt_residue_number = gt.loc[atom, 'residue_number']

        # Look for matching row in pred
        pred_row = pred.loc[(pred['atom_name'] == gt_atom_name) & (pred['residue_number'] == gt_residue_number) & (pred['residue_name'] == gt_residue_name)]
        if pred_row.empty != True:
            present_atoms.append(pred_row.index)

    # Select all rows in pred not present in gt
    total_atoms = list(pred.index)
    na_atoms_array = np.setdiff1d(total_atoms, present_atoms)
    na_atoms = sorted(na_atoms_array)

    # Create new pred data frame exclusively with atoms present in gt
    pred_trim = pred.drop(index=na_atoms)

    print(len(gt), len(pred_trim))
    # assert len(pred_trim) == len(gt)

    # print(f'Success! Creating trimmed files for {pdb}...')

    # new_gt = gt.to_csv(gt_out_path, sep='\t', index=False)
    # new_pred = pred_trim.to_csv(pred_out_path, sep='\t', index=False)
        # # pred_chain = pred.loc[pred['chain_id'] == gt_chain_id]
# # pred_residue = pred_chain.loc[pred_chain['residue_number'] == gt_residue_number]
# # pred_atom = pred_residue.loc[pred_residue['atom_name'] == gt_atom_name]
# pred_trim = pred.loc[(pred['chain_id'] == gt_chain_id) & (pred['residue_number'] == gt_residue_number) & (pred['atom_name'] == gt_atom_name)]