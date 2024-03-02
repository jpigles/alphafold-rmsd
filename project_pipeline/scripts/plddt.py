'''
A script for calculating average pLDDT scores of the entire protein and all annotated domains for given proteins.
'''

import pandas as pd
import main

# Define file paths
df = pd.read_csv(snakemake.input[0], sep='\t').astype('object')
af_auto = snakemake.input[1]
af_multi = snakemake.input[2]
af_single = snakemake.input[3]
cf_auto = snakemake.input[4]
cf_multi = snakemake.input[5]

# Calculate average pLDDT scores for AlphaFold2 models
af_auto_scores = main.calculate_pLDDT(df, af_auto)
af_multi_scores = main.calculate_pLDDT(df, af_multi)
af_single_scores = main.calculate_pLDDT(df, af_single)

# Calculate average pLDDT scores for ColabFold models
cf_auto_scores = main.calculate_pLDDT(df, cf_auto)
cf_multi_scores = main.calculate_pLDDT(df, cf_multi)

# Write results to file
af_auto_scores.to_csv(snakemake.output[0], sep='\t', index=False)
af_multi_scores.to_csv(snakemake.output[1], sep='\t', index=False)
af_single_scores.to_csv(snakemake.output[2], sep='\t', index=False)
cf_auto_scores.to_csv(snakemake.output[3], sep='\t', index=False)
cf_multi_scores.to_csv(snakemake.output[4], sep='\t', index=False)