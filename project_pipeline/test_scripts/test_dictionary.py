from pdbecif.mmcif_io import CifFileReader
from pdbecif.mmcif_tools import MMCIF2Dict
from pdbecif.mmcif_io import CifFileWriter
import pandas as pd
import numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch, Selection
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import test_utils
import requests
import shutil

'''
Merge the files.
'''
rmsds = pd.read_csv('./data/rmsds.tsv', sep='\t').astype('object')
classified_2 = pd.read_csv('./data/classified_files_2.tsv', sep='\t').astype('object')
pdbs = pd.read_csv('./data/proteins_by_pdb.tsv', sep='\t').astype('object')

rmsds = rmsds[['uniprot', 'pdb', 'complex_rmsd', 'percent_region_1', 'percent_region_2', '2_aligned', '2_comp']]
classified_2 = classified_2[['uniprot', 'pdb', 'conformation', 'state']]
pdbs = pdbs[['uniprot', 'pdb', 'region_1', 'region_2']]
print(len(rmsds))
print(len(classified_2))

classified_3 = pdbs.merge(rmsds, on=['uniprot', 'pdb'], how='right').drop_duplicates().reset_index(drop=True)
print(len(classified_3))
classified_3 = classified_3.merge(classified_2, on=['uniprot', 'pdb'], how='left').drop_duplicates().reset_index(drop=True)
print(len(classified_3))

'''
Retrieve the species from the cif file.
'''
for i in range(len(classified_3)):
    pdb = classified_3.loc[i, 'pdb']
    uniprot = classified_3.loc[i, 'uniprot']
    fp = './data/input/RCSB_cif/' + uniprot + '/' + pdb + '.cif'

    cfr = CifFileReader()
    cif_obj = cfr.read(fp, output='cif_wrapper')
    cif_data = list(cif_obj.values())[0]
    print(pdb)
    try:
        try:
            index = cif_data._struct_ref_seq.pdbx_db_accession.index(uniprot)
            org = cif_data._entity_src_gen.pdbx_gene_src_scientific_name[index]
        except:
            org = 'Manual Check Required'
        classified_3.loc[i, 'organism'] = org
    except AttributeError:
        try:
            index = cif_data._struct_ref_seq.pdbx_db_accession.index(uniprot)
            org = cif_data._entity_src_nat.pdbx_organism_scientific[index]
        except:
            org = 'Manual Check Required'
        classified_3.loc[i, 'organism'] = org

    date = cif_data._pdbx_database_status.recvd_initial_deposition_date
    classified_3.loc[i, 'date'] = str(date)


classified_3.to_csv('./data/classified_files_3.tsv', sep='\t', index=False)
# conf_cat.to_csv('./data/conf_cat.tsv', sep='\t', index=False)

# uniprot_cat.to_csv('./data/uniprot_cat.tsv', sep='\t', index=False)