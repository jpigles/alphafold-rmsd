#This will be used to scrape for all of the PDBs associated with a given Uniprot entry.
#But first, we have to get that list.

import csv

uniprot_id_list = []
dataset=[]

with open("combined_protein_list.csv") as protein_list:
    protein_list_dict = csv.DictReader(protein_list)
    for row in protein_list_dict:
        uniprot_id_list.append(row["Uniprot ID"])
        dataset.append(row)


uniprot_id_list_str = ""

for id in uniprot_id_list:
    uniprot_id_list_str = uniprot_id_list_str + id + " "

# print(uniprot_id_list_str)

import urllib.parse
import urllib.request

url = 'https://www.uniprot.org/uploadlists/'

params = {
'from': 'ACC+ID',
'to': 'PDB_ID',
'format': 'tab',
'query': uniprot_id_list_str
}

data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
   response = f.read()
pdbs_and_uniprot_str = response.decode('utf-8')

pdbs_and_uniprot_list = pdbs_and_uniprot_str.split("\n")
# print(pdbs_and_uniprot_list)
uniprot_only = pdbs_and_uniprot_list.pop(0)
# print(len(pdbs_and_uniprot_list))



#We want to go through every row in the csv file as a dictionary, check the value of the Uniprot ID key for that row against every item in the list I have made of the uniprot and PDB files

# Uncomment below to run the script
# for dictionary in dataset:
#     for item in pdbs_and_uniprot_list:
#         split_items = item.split("\t")
#         if split_items[0] == dictionary["Uniprot ID"]:
#                 dictionary["PDB"] = dictionary["PDB"] + " " + split_items[-1] 

# #This will create a new category, "Number of PDBs", that will gives us the total number of PDB files for a given entry.
# for dictionary in dataset:
#     dictionary["Number of PDBs"] = len(dictionary["PDB"].split())
# print(dataset)

# with open("combined_protein_list_pdb.csv", "w", newline="") as protein_list_pdb:
#     fields = ["Gene name", "IAS-range", "Accession Number", "Uniprot ID", "Referecne (PMID/uniprot)", "Protein Length", "Domain", "Domain-range", "Comments", "ENST", "Canonical_Source", "Protein name", 'Gene names', 'Protein length', "Species", "Autoinhibitory regions", "References", "Description (Refer to H and I)", "Evidence", "Structure?", "PDB", "Number of PDBs"]
#     output_writer = csv.DictWriter(protein_list_pdb, fieldnames=fields)

#     output_writer.writeheader()
#     for item in dataset:
#         output_writer.writerow(item)
