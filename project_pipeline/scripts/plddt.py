'''
A script for calculating average pLDDT scores of the entire protein and all annotated domains for given proteins.
'''

import pandas as pd
import main

# Define file paths
df_af_auto = pd.read_csv(snakemake.input[0], sep='\t').astype('object')
df_af_multi = pd.read_csv(snakemake.input[1], sep='\t').astype('object')
df_cf_auto = pd.read_csv(snakemake.input[2], sep='\t').astype('object')
df_cf_multi = pd.read_csv(snakemake.input[3], sep='\t').astype('object')
af_auto = snakemake.input[4]
af_multi = snakemake.input[5]
cf_auto = snakemake.input[6]
cf_multi = snakemake.input[7]

# Calculate average pLDDT scores for AlphaFold2 models
af_auto_scores = main.mean_plddt(df_af_auto, af_auto, unip_sub=False)
af_multi_scores = main.mean_plddt(df_af_multi, af_multi, unip_sub=False)

# Calculate average pLDDT scores for ColabFold models
cf_auto_scores = main.mean_plddt(df_cf_auto, cf_auto, unip_sub=True, fnt='cf_filename')
cf_multi_scores = main.mean_plddt(df_cf_multi, cf_multi, unip_sub=True, fnt='cf_filename')

# Write results to file
af_auto_scores.to_csv(snakemake.output[0], sep='\t', index=False)
af_multi_scores.to_csv(snakemake.output[1], sep='\t', index=False)
cf_auto_scores.to_csv(snakemake.output[2], sep='\t', index=False)
cf_multi_scores.to_csv(snakemake.output[3], sep='\t', index=False)