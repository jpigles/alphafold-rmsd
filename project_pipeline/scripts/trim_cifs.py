import pandas as pd
import main
import utils
import csv

gt_in_path = snakemake.input[0]
gt_trim_path = snakemake.input[1]
pred_in_path = snakemake.input[2]
pred_trim_path = snakemake.input[3]
best_cif_path = snakemake.input[4]
gt_file_format = snakemake.input[5]
pred_file_format = snakemake.input[6]

# Read in the reference dataframe
df_prot = pd.read_csv(snakemake.input[7], sep = '\t')

# Make the directories
utils.make_dirs(gt_in_path, gt_trim_path, pred_in_path, pred_trim_path, best_cif_path,
                 gt_format = gt_file_format, pred_format = pred_file_format)

# Trim the files to make sure any mutated residues are not accidentally counted in our domain completeness
trim_values, df_prot = main.trim_cifs(df_prot, gt_in_path, gt_trim_path, pred_in_path, pred_trim_path)

# Save the trim values
df_prot.to_csv(snakemake.output[0], sep = '\t', index = False)