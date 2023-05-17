'''
Analyze the quality of PDB files based on the percentage of residues in the inhibitory and active domains.
'''
import pandas as pd
import main


# Define the path of the CIF files and the reference dataframe
cif_in_path = 'data/input/RCSB_cif/'
df_prot = pd.read_csv(snakemake.input[0], sep = '\t')
cif_out_paths_list = ['data/input/RCSB_cif_80/', 'data/input/RCSB_cif_1_80/', 'data/input/RCSB_cif_2_80/', 'data/input/RCSB_cif_60/']

# Get the percentage of residues in the inhibitory and active domains
'''
Columns of data frame after this step: ['gene_name', 'uniprot', 'protein_length', 'region_1', 'region_2', 'region_1_len', 
                                 'region_2_len', 'pdb', 'pdb_length', 'resolution',
                                 'model', 'chain', 'auth_offset', 'pdb residues in region_1', 'pdb residues in region_2', 
                                 'percent_region_1', 'percent_region_2']
'''

print('Finding best files...')
df_prot = main.find_domain_completeness(df_prot, cif_in_path)

# Save several data frames with differentiating quality of domains: both domains with 80% of residues, domain 1 with 80% of residues,
# domain 2 with 80% of residues, and both domains with 60% of residues

df_list = main.save_domain_quality_files(df_prot, snakemake.output[0], snakemake.output[1], snakemake.output[2], snakemake.output[3], snakemake.output[4])

print('Found best files. Copying...')
for i in range(len(df_list)):
    copy_result = main.copy_best_files(df_list[i], cif_in_path, cif_out_paths_list[i])
    print(copy_result)

