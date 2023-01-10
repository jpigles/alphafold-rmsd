import requests

uniprot = 'Q16644'

url = f'https://www.ebi.ac.uk/pdbe/graph-api/uniprot/{uniprot}'

req = requests.get(url=url)

print(req.status_code)

req_json = req.json()

#To access our chain of interest, an example splice can be seen below

print(req_json[uniprot]['mappings'][0]['segments'][0]['chains'])