from biopandas.mmcif import PandasMmcif
import pandas as pd

df = pd.read_csv('./data/proteins_pdb_best.tsv', sep='\t').astype('object')

df_label = df.rename(columns={'Chain': 'Label_chain'})

auth_chains = []
for i in range(len(df)):

    # select PDB, file path, and chain
    pdb = df_label.loc[i, 'PDB ID']
    cif_path = (f'./data/input/RCSB_cif_best/{pdb}.cif')
    lbl_chain = df_label.loc[i, 'Label_chain']

    # Initialize PandasMmcif object and load mmcif file
    pmcif = PandasMmcif()
    _ = pmcif.read_mmcif(cif_path)
    pred = pmcif.df['ATOM']

    #Filter to only our label chain and select author
    auth_asym_df = pred.loc[pred['label_asym_id']==lbl_chain, 'auth_asym_id']
    auth_chain = auth_asym_df.iloc[0]
    auth_chains.append(auth_chain)

df_label.insert(12, 'Auth_chain', auth_chains)

df_label.to_csv('./data/proteins_pdb_best.tsv', sep='\t', index=False)