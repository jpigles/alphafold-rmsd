'''
A script for getting the RMSD between given AlphaFold2 models. In the context of this pipeline, one model is the public model 
and the other is the model from the AlphaFold2 pipeline.
'''

import os
import main
import pandas as pd

df = pd.read_csv('data/af2_structures.csv')

path1 = snakemake.input[0]
path2 = snakemake.input[1]
path3 = snakemake.input[2]

for index, row in df.iterrows():
    uniprot = row['uniprot']
    fn1 = f'F-{uniprot}-F1-model_v3.cif' # The public model
    fn2 = row['filename'] # The model from the AlphaFold2 pipeline
    complex_fn = fn2.split('-')[0] + '_comp.pdb' # eg P62826_U10-000_scores_rank_001_alphafold2_multimer_v2_model_1_seed_000.pdb -> P62826_U10_comp.pdb
    fp1 = os.path.join(path1, fn1)
    fp2 = os.path.join(path2, fn2)
    complex_out = os.path.join(path3, complex_fn)

    rmsd = main.complex_rmsd(fp1, fp2, complex_out)
    df.loc[index, 'rmsd'] = rmsd

df.to_csv(snakemake.output[0], index=False)