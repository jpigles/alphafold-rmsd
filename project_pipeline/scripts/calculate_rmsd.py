from pymol import cmd
import pandas as pd
import csv
import main

pdb_df = pd.read_csv(snakemake.input[0], sep='\t').astype('object')

rmsd_info = main.get_rmsds(pdb_df)

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