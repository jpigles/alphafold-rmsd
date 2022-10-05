# from Bio import SeqIO
# from functools import reduce
import os
from os.path import join

# n_g = 6
# fastas = []
# remove_x = False

# with open("./sample_data/input/idr/pdbs/1YCQ.pdb", 'r') as pdb_file: #Note that all of my files are .cif, so I will need to modify this. fn is .../inputidr_84/pdbs/{PDB_ID}.pdb.
#     for record in SeqIO.parse(pdb_file, 'pdb-atom'): #Here we change 'pdb-atom' to 'cif-atom'.
#         _ = record.id  #_ returns the value of the last executed expression value in the Python prompt/interpreter. So... what is that here? Well, now it's record.id
#         fastas.append(str(record.seq))
#         print(_)


# linker = 'G'*n_g #n_g equals 6 (Config[n_g])
# chain_start_resid_ids = {}
# linker_len = len(linker)

# seqs = fastas
# fasta = reduce(lambda acc, seq: acc + seq + linker, seqs, '')[:-n_g] #Sooooo, combine the sequences together, adding the linker to the end each time, then remove linker at the very end with the [:-n_g]

# acc, start_ids = 1, []
# for seq in seqs:
#     start_ids.append(acc) #for each file, add 1 as a start_id, then acc becomes 7 + len(seq)?
#     acc += len(seq) + n_g
# start_ids.append(len(fasta) + n_g + 1) # for ease of removal. Then we add the length of the combined chains + 7 to start_ids?
#         #I don't really understand what the point of this append is...

# chain_start_resid_ids['1YCQ'] = start_ids #dictionary of start_IDs associated with PDBs
# print(start_ids)
# print(chain_start_resid_ids)

# chain_start_ids = chain_start_resid_ids['1YCQ']

# n = len(chain_start_ids)
# for i in range(n-1): # for each chain
#     resid_lo = chain_start_ids[i]
#     # resid_hi is residue id of first 'g' in pred
#     resid_hi = chain_start_ids[i+1] - linker_len
#     # print(resid_lo)
#     # print(resid_hi)
pdb = 'IE3C'
snakemake = 'project_pipeline/sample_data/output/idr_af_full/poly_g_6'
output_dir = os.mkdir(snakemake + f'/{pdb}.fasta')
output_path = join(snakemake, f'{pdb}.fasta', 'ranked_0.pdb')
stream = os.popen(f'python ./env/lib/python3.9/site-packages/pdbtools/pdb_fromcif.py project_pipeline/sample_data/input/idr/cif/F-B8XX90-F1-model_v3.cif > {output_path}')

# path = join('data/input/poly_g_6', 'pdb.fasta', 'uniprot.pdb')
# print(path)