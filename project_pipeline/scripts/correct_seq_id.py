from biopandas.pdb import PandasPdb
import pandas as pd
from pdbecif.mmcif_io import CifFileReader

# Load list of our best pdb files
pdbs_df = pd.read_csv('./data/proteins_pdb_no_offset.tsv', sep = '\t').astype('object')

# Functions start here
def get_offset(fp, pdb):
    # initiate reader object
    cfr = CifFileReader()
    cif_obj = cfr.read(fp, output='cif_wrapper')
    cif_data = list(cif_obj.values())[0]
    
    # Extract the auth_seq start and db_seq start (from Uniprot) from the cif file
    auth_start = int(cif_data._struct_ref_seq.pdbx_auth_seq_align_beg[0])
    unp_start = int(cif_data._struct_ref_seq.db_align_beg[0])
    offset = auth_start - unp_start

    print(f'Offset for {pdb}: {offset}')
    return offset

def fix_seq_id(pdb, in_fp, out_fp, chain, offset):
    # initiate Pandas object
    ppdb = PandasPdb()
    _ = ppdb.read_pdb(in_fp)
    pred = ppdb.df['ATOM']

    # Replace residue numbers in our chain of interest
    if offset == 0:
        return f'No fix needed for {pdb}'
    else:
        for i in range(len(pred)):
            if pred.loc[i, 'chain_id']==chain:
                res_num = pred.loc[i, 'residue_number']
                new_res_num = res_num - offset
                pred.loc[i, 'residue_number'] = new_res_num
            else:
                continue
    
    ppdb.to_pdb(path=out_fp, records=None, gz=False, append_newline=True)

    return f'Successfully fixed {pdb}'

# Get offset between author and uniprot seq id
offsets = []

# Run through main list and fix pdb files
for i in range(len(pdbs_df)):
    # Define pdbs to work with and auth_chains in pdb files
    pdb_id = pdbs_df.loc[i, 'PDB ID']
    auth_chain = pdbs_df.loc[i, 'Auth_chain']

    # Designate file locations
    cif_path = f'./data/input/RCSB_cif_best/{pdb_id}.cif'
    pdb_path = f'./data/input/RCSB/offset_pdbs/{pdb_id}.pdb'
    out_path = f'./data/input/RCSB/pdbs/{pdb_id}.pdb'

    # Get offset
    offset= get_offset(cif_path, pdb_id)
    offsets.append(offset)

    # fix pdb sequence ids
    fixed_pdb = fix_seq_id(pdb_id, pdb_path, out_path, auth_chain, offset)
    print(fixed_pdb)

# Add offsets to file
pdbs_df.insert(18, 'Auth_offset', offsets)
pdbs_df.to_csv('./data/proteins_pdb_best.tsv', sep = '\t', index=False)