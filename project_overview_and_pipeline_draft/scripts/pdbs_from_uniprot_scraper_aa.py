#2022-08-23: This is a modification of the original scraper. However, it doesn't seem to work because of a URL error.

import csv

# uniprot_id_list = []
# dataset=[]

# with open("/home/bjechow/Documents/gsponer_lab/autoinhibition_protein_data/project_overview_and_pipeline_draft/data/proteins.tsv", newline="") as protein_list:
#     protein_list_dict = csv.DictReader(protein_list, delimiter="\t")
#     for row in protein_list_dict:
#         uniprot_id_list.append(row["Uniprot_ID"])
#         dataset.append(row)


# uniprot_id_list_str = ""

# for id in uniprot_id_list:
#     uniprot_id_list_str = uniprot_id_list_str + id + " "

# # print(uniprot_id_list_str)

# import urllib.parse
# import urllib.request

# url = 'https://www.uniprot.org/uploadlists/'

# params = {
# 'from': 'ACC+ID',
# 'to': 'PDB_ID',
# 'format': 'tab',
# 'query': uniprot_id_list_str
# }

# data = urllib.parse.urlencode(params)
# data = data.encode('utf-8')
# req = urllib.request.Request(url, data)
# with urllib.request.urlopen(req) as f:
#    response = f.read()
# pdbs_and_uniprot_str = response.decode('utf-8')

# pdbs_and_uniprot_list = pdbs_and_uniprot_str.split("\n")
# # print(pdbs_and_uniprot_list)
# uniprot_only = pdbs_and_uniprot_list.pop(0)
# print(len(pdbs_and_uniprot_list))



# #We want to go through every row in the csv file as a dictionary, check the value of the Uniprot ID key for that row against every item in the list I have made of the uniprot and PDB files

# # Uncomment below to run the script
# for dictionary in dataset:
#     for item in pdbs_and_uniprot_list:
#         split_items = item.split("\t")
#         if split_items[0] == dictionary["Uniprot_ID"]:
#                 dictionary["PDB"] = dictionary["PDB"] + "," + split_items[-1] 

# #This will create a new category, "Number of PDBs", that will give us the total number of PDB files for a given entry.
# for dictionary in dataset:
#     dictionary["Number_of_PDBs"] = len(dictionary["PDB"].split())
# # print(dataset)

# with open("/home/bjechow/Documents/gsponer_lab/autoinhibition_protein_data/project_overview_and_pipeline_draft/data/proteins_pdb.tsv", "w", newline="") as protein_list_pdb:
#     fields = ["Gene_name", "Uniprot_ID", "Protein_length", "region_1", "region_2", "PDB", "Number_of_PDBs"]
#     output_writer = csv.DictWriter(protein_list_pdb, fieldnames=fields)

#     output_writer.writeheader()
#     for item in dataset:
#         output_writer.writerow(item)


# The below section is for clearing out white space at the start of each PDB file list in our csv.

pdb_id_list = []
dataset=[]
trimmed_pdb_id_list = []

with open("/home/bjechow/Documents/gsponer_lab/autoinhibition_protein_data/project_overview_and_pipeline_draft/data/proteins.tsv", newline="") as protein_list:
    protein_list_dict = csv.DictReader(protein_list, delimiter="\t")
    for row in protein_list_dict:
        pdb_id_list.append(row["PDB"])
        dataset.append(row)

for dictionary in dataset:
    dictionary["PDB"] = dictionary["PDB"].strip()
    # print(dictionary["PDB"])

with open("/home/bjechow/Documents/gsponer_lab/autoinhibition_protein_data/project_overview_and_pipeline_draft/data/proteins.tsv", "w", newline="") as protein_list_pdb:
    fields = ["Gene_name", "Uniprot_ID", "Protein_length", "region_1", "region_2", "PDB"]
    output_writer = csv.DictWriter(protein_list_pdb, fieldnames=fields)

    output_writer.writeheader()
    for item in dataset:
        output_writer.writerow(item)
