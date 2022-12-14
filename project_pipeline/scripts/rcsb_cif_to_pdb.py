import pandas as pd
import os

from os.path import join

#Create the dataframe
prot_df = pd.read_csv(snakemake.input[0], sep='\t').astype('object')

# Set up the appropriate pdb/uniprot relationship
for i in range(len(prot_df)):
    pdb = prot_df.loc[i, 'PDB ID']
    uniprot = prot_df.loc[i, 'Uniprot_ID']
    
    #Make output dir
    try:
        output_dir = os.mkdir(snakemake.output[0])
    except:
        pass
    
    #Define paths
    input_path = join(snakemake.input[1], f'{pdb}.cif')
    output_path = join(snakemake.output[0], f'{pdb}.pdb')

    # Convert the files
    print(f'Converting {pdb}.cif to .pdb...')
    stream = os.popen(f'python ../env/lib/python3.9/site-packages/pdbtools/pdb_fromcif.py {input_path} > {output_path}')