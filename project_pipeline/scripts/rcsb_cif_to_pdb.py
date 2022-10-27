import pandas as pd
import os

from os.path import join

#Create the dataframe
af_df = pd.read_csv(snakemake.input[0], sep='\t').astype('object')

# Set up the appropriate pdb/uniprot relationship
for i in range(len(af_df)):
    pdb = af_df.loc[i, 'PDB ID']
    uniprot = af_df.loc[i, 'Uniprot_ID']
    
    #Make output dir
    output_dir = os.mkdir(snakemake.output[0])

    #Define paths
    input_path = join(snakemake.input[1], f'{pdb}.cif')
    output_path = join(snakemake.output[0], f'{pdb}.pdb')

    #Convert the files
    stream = os.popen(f'python ./env/lib/python3.9/site-packages/pdbtools/pdb_fromcif.py {input_path} > {output_path}')