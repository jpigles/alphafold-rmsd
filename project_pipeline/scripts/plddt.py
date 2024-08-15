'''
A script for calculating average pLDDT scores of the entire protein and all annotated domains for given proteins.
'''

import pandas as pd
import main

# Define file paths
df_ai = pd.read_csv(snakemake.input[0], sep='\t').astype('object')
df_md = pd.read_csv(snakemake.input[1], sep='\t').astype('object')  # TODO: Pass mean_plddt one dataframe containing only uniprots and af_filenames,
df_sd = pd.read_csv(snakemake.input[2], sep='\t').astype('object')   # and another (unchanged) dataframe containing cf_filenames. 
ai_fd_fp = snakemake.input[3]
md_fd_fp = snakemake.input[4]
sd_fd_fp = snakemake.input[5]
ai_c_fp = snakemake.input[6]
md_c_fp = snakemake.input[7]

# Collapse autoinhibitory and multi-domain dataframes by dropping clusters and pdb info
df_ai_fd = df_ai[['uniprot', 'region_1', 'region_2', 'af_filename', 'cf_filename']].drop_duplicates().reset_index(drop=True)
df_md_fd = df_md[['uniprot', 'region_1', 'region_2', 'af_filename', 'cf_filename']].drop_duplicates().reset_index(drop=True)

# Calculate average pLDDT scores for AlphaFold2 models
ai_fd_scores = main.mean_plddt(df_ai_fd, ai_fd_fp)
md_fd_scores = main.mean_plddt(df_md_fd, md_fd_fp)
sd_fd_scores = main.mean_plddt_single_domain(df_sd, sd_fd_fp)

# Calculate average pLDDT scores for ColabFold models
ai_c_scores = main.mean_plddt(df_ai, ai_c_fp, fnt='cf_filename')
md_c_scores = main.mean_plddt(df_md_fd, md_c_fp, fnt='cf_filename')

# Write results to file
ai_fd_scores.to_csv(snakemake.output[0], sep='\t', index=False)
md_fd_scores.to_csv(snakemake.output[1], sep='\t', index=False)
sd_fd_scores.to_csv(snakemake.output[2], sep='\t', index=False)
ai_c_scores.to_csv(snakemake.output[3], sep='\t', index=False)
md_c_scores.to_csv(snakemake.output[4], sep='\t', index=False)