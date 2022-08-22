import csv

#Contains alphafold uniprots
alphafold_uniprot_list = []

#Contains curated uniprots
curated_uniprot_list = []

#Where to store the present ids
present_ids = []

#Where to store the missing ids
missing_ids = []

#Store alphafold uniprot IDs in a dictionary
with open(r"/home/bjechow/Documents/gsponer_lab/autoinhibition_protein_data/curated_protein_lists/determine_missing_uniprot/autoinhibited_proteins_uniprot_from_alphafold.csv", newline="") as alphafold:
    alphafold_dict = csv.DictReader(alphafold)
    for row in alphafold_dict:
        alphafold_uniprot_list.append(row)
        if row["Uniprot"] not in present_ids:
            present_ids.append(row["Uniprot"])

#Store curated uniprot IDs in a dictionary
with open(r"/home/bjechow/Documents/gsponer_lab/autoinhibition_protein_data/curated_protein_lists/determine_missing_uniprot/autoinhibited_proteins_uniprot.csv", newline="") as curated:
    curated_dict = csv.DictReader(curated)
    for row in curated_dict:
        curated_uniprot_list.append(row)
        if row["Uniprot"] not in present_ids:
            missing_ids.append(row)

# Uncomment below to run the script.
# with open("./curated_protein_lists/determine_missing_uniprot/missing_uniprots_in_alphafold.csv", "w", newline="") as missing_uniprots:
#     fields = ["Uniprot"]
#     output_writer = csv.DictWriter(missing_uniprots, fieldnames=fields)

#     output_writer.writeheader()
#     for item in missing_ids:
#         output_writer.writerow(item)