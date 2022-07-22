import csv

#Contains all proteins
master_uniprot_list = []
#Contains all repeated Uniprot IDs
repeats_list = []
no_uniprot_repeats = []
#contains all information from Brooks file
master_brooks_list = []
#contains all information from Jorge file
master_jorge_file = []

combined_protein_file = []


# This will import all of the Uniprot ID's from the Autoinhibited Proteins file
with open(r"C:\Users\bkadu\Documents\Python_Projects\Compiling_Protein_Lists\Final_compilation\Autoinhibited proteins (more added by Brooks).csv", newline="") as brooks_file:
    brooks_dict = csv.DictReader(brooks_file)
    for row in brooks_dict:
        master_brooks_list.append(row)
        if row["Uniprot ID"] not in master_uniprot_list:
            master_uniprot_list.append(row["Uniprot ID"])

#This will import all of the Uniprot ID's from Jorge's file
with open(r"C:\Users\bkadu\Documents\Python_Projects\Compiling_Protein_Lists\Final_compilation\Jorge_Algonquin_Autoinhibition_Proteins_w_Structure.csv", newline="") as jorge_file:
    jorge_dict = csv.DictReader(jorge_file)
    for row in jorge_dict:
        master_jorge_file.append(row)
        if row["Uniprot ID"] not in master_uniprot_list:
            master_uniprot_list.append(row["Uniprot ID"])
        else:
            repeats_list.append(row["Uniprot ID"])

for bdict in master_brooks_list:
    for jdict in master_jorge_file:
        if bdict["Uniprot ID"] == jdict["Uniprot ID"]:
            combined_protein = bdict | jdict
            combined_protein_file.append(combined_protein)


for id in master_uniprot_list:
    if id not in repeats_list:
        no_uniprot_repeats.append(id)


for bdict in master_brooks_list:
    if bdict["Uniprot ID"] in no_uniprot_repeats:
        combined_protein_file.append(bdict)

for jdict in master_jorge_file:
    if jdict["Uniprot ID"] in no_uniprot_repeats:
        combined_protein_file.append(jdict)

print(len(repeats_list))

# Uncomment this below if you want to run the script again.
with open("combined_protein_list.csv", "w", newline="") as combined_list:
    fields = ["Gene name", "IAS-range", "Accession Number", "Uniprot ID", "Referecne (PMID/uniprot)", "Protein Length", "Domain", "Domain-range", "Comments", "ENST", "Canonical_Source", "Protein name", "Gene names", "Protein length", "Species", "Autoinhibitory regions", "References", "Description (Refer to H and I)", "Evidence", "Structure?", "PDB", '']
    output_writer = csv.DictWriter(combined_list, fieldnames=fields)

    output_writer.writeheader()
    for item in combined_protein_file:
        output_writer.writerow(item)