'''
This is a script for taking the trimmed PDB files and splitting them into two chains
for use in DockQ. The autoinhibitory region will become chain B, while the functional domain
and the rest of the protein will remain chain A.
'''
import pandas as pd
import utils
import main

df = pd.read_csv(snakemake.input[0], sep='\t')
gt_path_in = snakemake.input[1]
pred_path_in = snakemake.input[2]
gt_path_out = snakemake.output[0]
pred_path_out = snakemake.output[1]

# Make directories
utils.make_dirs(gt_path_in, pred_path_in, gt_path_out, pred_path_out)

# Split chains
main.split_chains(df, gt_path_in, pred_path_in, gt_path_out, pred_path_out)