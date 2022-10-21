# from Bio import SeqIO
# from functools import reduce
import requests
import os
from os.path import join
import pandas
from scripts.mutation_enrichment import string2range

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
# pdb = 'IE3C'
# snakemake = 'project_pipeline/sample_data/output/idr_af_full/poly_g_6'
# output_dir = os.mkdir(snakemake + f'/{pdb}.fasta')
# output_path = join(snakemake, f'{pdb}.fasta', 'ranked_0.pdb')
# stream = os.popen(f'python ./env/lib/python3.9/site-packages/pdbtools/pdb_fromcif.py project_pipeline/sample_data/input/idr/cif/F-B8XX90-F1-model_v3.cif > {output_path}')

# path = join('data/input/poly_g_6', 'pdb.fasta', 'uniprot.pdb')
# print(path)

# df_prot = pandas.read_csv('sample_data/proteins_pdb_sample.csv', sep = ',').astype('object')


# aa_range = '166-309,482-618'
# def string2range(x):

#     if ',' in x:
#             list_temp = x.split(sep = ',') #list_temp = ['123-456,' '789-1111']
#             # print(f'This is list_temp: {list_temp}')
#             for y in range(len(list_temp)): 
#                 list_temp[y] = list_temp[y].split(sep = '-') #list_temp[y] = [['123', '456'], ['789', '1111']]
#             for y in range(len(list_temp)): 
#                 for x in range(len(list_temp[y])):
#                     list_temp[y][x] = int(list_temp[y][x]) #turns each list item into an integer
#                     # print(f'This is list_temp[y][x]: {list_temp[y][x]}')

#             # Make a range object with the bounds of the range. Note to the 
#             # end a 1 has to be added in order to include the last position in the range
#             for y in range(len(list_temp)): #[1, 2] where 1=[123, 456] and 2=[789, 1111]
#                 for x in range(len(list_temp[y])): #[123, 456]       
#                     list_temp[y] = list(range(list_temp[y][x], list_temp[y][x+1]+1)) #list_temp[0][0] = [123], list_temp[0][0+1]+1 or [456] + 1 = [457]
#                     # print(f'This is list_temp[y]:{list_temp[y]}')
#                     break

#             # print(list(set([item for sublist in list_temp for item in sublist])))
#             return list(set([item for sublist in list_temp for item in sublist]))

#         # Handle instances with only one range
#     else:
#         list_temp = x.split(sep = '-')
#         for y in range(len(list_temp)):
#             list_temp[y] = int(list_temp[y]) #

#         # Make a range object with the bounds of the region. Note to the 
#         # end a 1 has to be added in order to include the last position in the range
#         return list(range(list_temp[0], list_temp[1]+1))



# df_prot['region_1 search'] = df_prot['region_1'].apply(lambda x: string2range(x))
# df_prot['region_2 search'] = df_prot['region_2'].apply(lambda x: string2range(x)) 


# print(df_prot['region_2 search'])

# def prune_extra_chains(pdb_ids):
#     pdb_ids_only = []
#     single_chain_ids = set()
#     for pdb_id in pdb_ids:
#         id_only = pdb_id[:4]
#         lowercase_id = id_only.lower()
#         pdb_ids_only.append(lowercase_id)
#     for pdb_id in pdb_ids_only:
#         if pdb_ids_only.count(pdb_id) != 1:
#             pdb_ids_only.remove(pdb_id)
#         else:
#             single_chain_ids.add(pdb_id)
#     return single_chain_ids

# id_str = '4Y07.A 5TJ7.A 5TJ7.B 5TJ7.C 5TJ7.D 5TJ8.A 5TJQ.A 6J1Z.A 6RSS.A '

# id_list = id_str.strip().split(sep=' ')

# pruned_list = prune_extra_chains(id_list)

# print(pruned_list)