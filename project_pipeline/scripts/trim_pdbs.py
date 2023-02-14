from biopandas.pdb import PandasPdb
import pandas as pd
import numpy as np
import os
import csv

pdb_list = pd.read_csv('./data/proteins_pdb_best.tsv', sep='\t').astype('object')

too_short_trims = []
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
    gt_path = f'./data/input/RCSB/pdbs/{pdb}.pdb'
    pred_path = f'./data/output/RCSB_af_full/poly_g_6_fasta/{pdb}.fasta/ranked_0.pdb'
    gt_out_path = f'./data/input/RCSB/pdbs_trim/{pdb}.pdb'
    pred_out_path = f'./data/output/RCSB_af_full/af_trim/{pdb}.fasta/ranked_0.pdb'

    print(f'Trying {pdb}...')

    # gt model
    ppdb = PandasPdb()
    gt_file = ppdb.read_pdb(gt_path)
    gt_all_chains = gt_file.df['ATOM']
    gt = gt_all_chains[gt_all_chains['chain_id'] == auth_chain]
    gt = gt.reset_index()

    pred_file = ppdb.read_pdb(pred_path)
    pred = pred_file.df['ATOM']

    print(len(gt), len(pred))

    # Create list of rows that are present in both gt and pred based on row indices of pred
    present_atoms_pred = []
    extra_atoms_gt = []

    # Define atom_names for hydrogens
    hydrogens = ['HA', 'HB1', 'HB2', 'HB3', 'H', 'HA2', 'HA3', 'HG2', 'HG3', 'HD2', 'HD3', 'HE1', 'HE2',
                'HE3', 'HB', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23', 'HE', 'HH11', 'HH12', 'HH21',
                'HH22', 'HE21', 'HE22', 'HD1', 'HZ', 'HH', 'HG1', 'HD21', 'HD22', 'HG', 'HD11', 'HD12', 
                'HD13', 'HD23', 'HZ1', 'HZ2', 'HZ3', 'HH2']
    # Define possible alternate conformations
    alt_locations = ['B', 'C', 'D', 'E']

    for atom in range(len(gt)):
        # Define rows to be skipped (hydrogens or alternate conformations)
        if gt.loc[atom, 'atom_name'] in hydrogens or gt.loc[atom, 'alt_loc'] in alt_locations:
            extra_atoms_gt.append(atom)
            continue

        # Define minimum  parameters to select unique rows
        gt_atom_name = gt.loc[atom, 'atom_name']
        gt_residue_name = gt.loc[atom, 'residue_name']
        gt_residue_number = gt.loc[atom, 'residue_number']

        # Look for matching row in pred
        pred_row = pred.loc[(pred['atom_name'] == gt_atom_name) & (pred['residue_number'] == gt_residue_number) & (pred['residue_name'] == gt_residue_name)]
        if pred_row.empty != True:
            present_atoms_pred.append(pred_row.index)
        else:
            extra_atoms_gt.append(atom)

    # Select all rows in pred not present in gt
    total_atoms = list(pred.index)
    na_atoms_array = np.setdiff1d(total_atoms, present_atoms_pred)
    na_atoms = sorted(na_atoms_array)

    # Create new pred data frame exclusively with atoms present in gt
    pred_trim = pred.drop(index=na_atoms)
    gt_trim = gt.drop(index=extra_atoms_gt)

    print('Length of gt: ' + str(len(gt_trim)) + ', Length of pred: ' + str(len(pred_trim)))

    try:
        assert len(pred_trim) == len(gt_trim)
    except AssertionError:
        gt_sim = gt_trim.drop(['atom_number', 'x_coord', 'y_coord', 'z_coord', 'occupancy', 'b_factor', 'segment_id', 'element_symbol', 'charge', 'line_idx'], axis=1)
        pred_sim = pred_trim.drop(['atom_number', 'x_coord', 'y_coord', 'z_coord', 'occupancy', 'b_factor', 'segment_id', 'element_symbol', 'charge', 'line_idx'], axis=1)
        diff = pd.concat([gt_sim, pred_sim]).drop_duplicates(keep=False)
        diff.to_csv('./data/AssertionError.tsv', sep='\t')
        print(diff)
        print('AssertionError! Check file')
        break

    if len(gt_trim) < 1500:
        too_short_trims.append(pdb)

    print(f'Success! Creating trimmed files for {pdb}...')
    # Make directory for specific pdb
    try:
        new_gt = gt_trim.to_csv(gt_out_path, sep='\t', index=False)
        new_pred = pred_trim.to_csv(pred_out_path, sep='\t', index=False)
    except OSError:
        af_dir = os.mkdir(f'./data/output/RCSB_af_full/af_trim/{pdb}.fasta/')
        new_pred = pred_trim.to_csv(pred_out_path, sep='\t', index=False)

with open('./data/wrong_offsets.tsv', 'w') as file:
    fields = ['Incorrect_pdbs']
    writer = csv.writer(file)
    
    writer.writerow(fields)
    writer.writerows(too_short_trims)