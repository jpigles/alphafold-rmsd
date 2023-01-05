from biopandas.mmcif import PandasMmcif
import pandas as pd

df = pd.read_csv('../sample_data/proteins_pdb_best.csv', sep=',').astype('object')

df_label = df.rename(columns={'Chain': 'Label_chain'})

auth_chains = []
for i in range(len(df)):

    # select PDB, file path, and chain
    pdb = df_label.loc[i, 'PDB']
    cif_path = (f'../data/input/RCSB_cif/{pdb}.cif')
    lbl_chain = df_label.loc[i, 'Label_chain']

    # Initialize PandasMmcif object and load mmcif file
    pmcif = PandasMmcif()
    pred = pmcif.read_mmcif(cif_path)

    #Filter to only our label chain and select author
    df_one_chain = df_label['label_asym_id']==lbl_chain
    auth_chain = df_one_chain.loc[0, 'auth_asym_id']
    auth_chains = auth_chains.append(auth_chain)

df_auth = df_label.insert(13, 'Auth_chain', auth_chains)

df_auth.to_csv('../sample_data/proteins_auth_chains.tsv', sep='\t', index=False)