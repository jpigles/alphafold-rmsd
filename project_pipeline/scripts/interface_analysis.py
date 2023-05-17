'''
Analyze the quality of PDB files based on the percentage of residues in the inhibitory and active domains.
'''
import pandas as pd
import main


# Define the path of the CIF files and the reference dataframe
cif_path = 'data/input/RCSB_cif/'
df_prot = pd.read_csv(snakemake.input[0], sep = '\t')

# Get the percentage of residues in the inhibitory and active domains
'''
Columns of data frame after this step: ['gene_name', 'uniprot', 'protein_length', 'region_1', 'region_2', 'region_1_len', 
                                 'region_2_len', 'pdb', 'pdb_length', 'resolution',
                                 'model', 'chain', 'auth_offset', 'pdb residues in region_1', 'pdb residues in region_2', 
                                 'percent_region_1', 'percent_region_2']
'''
df_prot = main.find_domain_completeness(df_prot, cif_path)

# Save several data frames with differentiating quality of domains: both domains with 80% of residues, domain 1 with 80% of residues,
# domain 2 with 80% of residues, and both domains with 60% of residues

main.save_domain_quality_files(df_prot, snakemake.output[0], snakemake.output[1], snakemake.output[2], snakemake.output[3], snakemake.output[4])