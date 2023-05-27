from pymol import cmd
import pandas as pd
import csv
import main

gt_in_path = './data/input/RCSB_cif_best/'
gt_trim_path = './data/input/RCSB_cif_trim/'
pred_in_path = './data/input/Alphafold_cif/'
pred_trim_path = './data/input/Alphafold_cif_trim/'
complex_path = './data/output/complexes/'
df = pd.read_csv(snakemake.input[0], sep='\t').astype('object')

# First we trim the files
trim_values = main.trim_cifs(df, gt_in_path, gt_trim_path, pred_in_path, pred_trim_path)

# Save the trim values
with open('./data/trim_values.tsv', 'w') as file:
    fields = ['pdb', 'gt_len', 'gt_trim_len', 'pred_len', 'pred_trim_len', 'gt_perc', 'trim_perc']
    writer = csv.DictWriter(file, fieldnames=fields, delimiter='\t')
    
    writer.writeheader()
    for item in trim_values:
        writer.writerow(item)

# Calculate and save the rmsd info
rmsd_info = main.get_rmsds(df, gt_trim_path, pred_trim_path, complex_path)

with open(snakemake.output[0], 'w') as file:
    fields = ['UniProt', 'PDB', 'complex_rmsd', '1.0_aligned', '1.0_comp',
                '1.1_aligned', '1.1_comp', '1.2_aligned', '1.2_comp', '2.0_aligned', '2.0_comp',
                '2.1_aligned', '2.1_comp', '2.2_aligned', '2.2_comp', '2.3_aligned', '2.3_comp',
                'Percent residues in region_1', 'Percent residues in region_2', '1_aligned', '1_comp',
                '2_aligned', '2_comp']
    writer = csv.DictWriter(file, fieldnames=fields, delimiter='\t')
    
    writer.writeheader()
    for item in rmsd_info:
        writer.writerow(item)