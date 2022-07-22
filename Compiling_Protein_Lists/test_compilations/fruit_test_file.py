import csv

combined_fruit_name_list = []
repeated_names = []
no_repeated_names = []


fruit2_list = []
fruit1_list = []
combined_fruit_list = []

with open("fruit_2.csv") as fruit2:
    fruit2_dict = csv.DictReader(fruit2)
    for row in fruit2_dict:
        fruit2_list.append(row)
        if row["ï»¿Fruit"] not in combined_fruit_name_list:
            combined_fruit_name_list.append(row["ï»¿Fruit"])


# print(fruit2_list)

with open("fruits_1.csv") as fruit1:
    fruit1_dict = csv.DictReader(fruit1)
    for row in fruit1_dict:
        fruit1_list.append(row)
        if row["ï»¿Fruit"] not in combined_fruit_name_list:
            combined_fruit_name_list.append(row["ï»¿Fruit"])
        else:
            repeated_names.append(row["ï»¿Fruit"])

for dict2 in fruit2_list:
    for dict1 in fruit1_list:
        if dict1["ï»¿Fruit"] == dict2["ï»¿Fruit"]:
            combined_fruit = dict1 | dict2
            combined_fruit_list.append(combined_fruit)

# for dict2 in fruit2_list:
#     if dict2["ï»¿Fruit"] in combined_fruit_name_list

for name in combined_fruit_name_list:
    if name not in repeated_names:
        no_repeated_names.append(name)


for dict2 in fruit2_list:
    if dict2["ï»¿Fruit"] in no_repeated_names:
        combined_fruit_list.append(dict2)

for dict1 in fruit1_list:
    if dict1["ï»¿Fruit"] in no_repeated_names:
        combined_fruit_list.append(dict1)


with open("combined_fruit.csv", "w") as combined_fruit:
    fields = ["ï»¿Fruit", "Price", "Quality", "Source", "Age"]
    output_writer = csv.DictWriter(combined_fruit, fieldnames=fields)

    output_writer.writeheader()
    for item in combined_fruit_list:
        output_writer.writerow(item)