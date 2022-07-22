import csv

test_list = ['From\tTo', 'Q9UM73\t2KUP', 'Q9UM73\t2KUQ', 'Q9UM73\t2XB7', 'Q9UM73\t2XBA', 'Q9UM73\t2XP2', 'Q9UM73\t2YFX', 'Q9UM73\t2YHV', 'Q9UM73\t2YJR', 'Q9UM73\t2YJS', 'Q9UM73\t2YS5']
uniprot_only = test_list.pop(0)


#We want to go through every row in the csv file as a dictionary, check the value of the Uniprot ID key for that row against every item in the list I have made of the uniprot and PDB files
dataset = []

with open("small_combined_protein_dataset.csv",) as small_dataset:
    small_dataset_dict = csv.DictReader(small_dataset, skipinitialspace=True)
    for row in small_dataset_dict:
        dataset.append(row)

# print(dataset)

for dictionary in dataset:
    for item in test_list:
        split_items = item.split("\t")
        if split_items[0] == dictionary["Uniprot ID"]:
                dictionary["PDB"] = dictionary["PDB"] + split_items[-1] + " "

print(dataset)

with open("small_combined_protein_dataset_add_pdb.csv", "w", newline="") as small_dataset:
    fields = ["Gene name", "IAS-range", "Accession Number", "Uniprot ID", "Referecne (PMID/uniprot)", "Protein Length", "Domain", "Domain-range", "Comments", "ENST", "Canonical_Source", "Protein name", 'Gene names', 'Protein length', "Species", "Autoinhibitory regions", "References", "Description (Refer to H and I)", "Evidence", "Structure?", "PDB"]
    output_writer = csv.DictWriter(small_dataset, fieldnames=fields)

    output_writer.writeheader()
    for item in dataset:
        output_writer.writerow(item)


    
    



# with open("small_combined_protein_dataset_pdb.csv", "w") as expand_dataset:
#     fields = ["Gene name", "IAS-range", "Accession Number", "Uniprot ID", 'Referecne (PMID/uniprot)', 'Protein Length', 'Domain', 'Domain-range', 'Comments', 'ENST', 'Canonical_Source', 'Protein name', 'Gene names', 'Protein length', 'Species', 'Autoinhibitory regions', 'References', 'Description (Refer to H and I)', 'Evidence', 'Structure?', 'PDB']
#     expand_dataset_dict = csv.DictWriter(expand_dataset, fieldnames=fields)
    


#open my CSV file in append mode
#Take the items in test list and split them
#for row in combined protein dictionary:
#   for item in test_list:
#       split_items = item.split("\t")
#       if split_items[0] == row["Uniprot ID"]:
#           row["PDB"] = row["PDB"] + split_items[-1] + " "
