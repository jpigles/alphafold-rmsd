import pandas as pd
import main
import utils
import csv

gt_in_path = 'data/input/RCSB_cif/'
gt_trim_path = './data/input/RCSB_cif_trim/'
pred_in_path = './data/input/Alphafold_cif/'
pred_trim_path = './data/input/Alphafold_cif_trim/'
best_cif_path = 'data/input/RCSB_cif_best/'

# Read in the reference dataframe
df_prot = pd.read_csv(snakemake.input[0], sep = '\t')

# Make the directories
utils.make_dirs([gt_in_path, gt_trim_path, pred_in_path, pred_trim_path, best_cif_path])

# Trim the files to make sure any mutated residues are not accidentally counted in our domain completeness
trim_values, df_prot = main.trim_cifs(df_prot, gt_in_path, gt_trim_path, pred_in_path, pred_trim_path)

# Save the trim values
df_prot.to_csv(snakemake.output[0], sep = '\t', index = False)