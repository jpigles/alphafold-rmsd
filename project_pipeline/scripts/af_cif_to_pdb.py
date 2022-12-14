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
    # print(snakemake.output[0] + f'/{pdb}.fasta/')
    if not os.path.exists('/data/output/RCSB_af_full/poly_g_6' + '/' + pdb + '.fasta/'):
        try:
            original_umask = os.umask(0)
            os.makedirs('data/output/RCSB_af_full/poly_g_6' + '/' + pdb + '.fasta/', 0o0770)
        except:
            pass
        finally:
            os.umask(original_umask)
    else:
        pass

    #Define paths
    input_path = join(snakemake.input[1], f'F-{uniprot}-F1-model_v3.cif')
    output_path = join(snakemake.output[0], f'{pdb}.fasta', 'ranked_0.pdb')

    #Convert the files
    print(f'Converting {uniprot}.cif to .pdb...')
    stream = os.popen(f'python ../env/lib/python3.9/site-packages/pdbtools/pdb_fromcif.py {input_path} > {output_path}')