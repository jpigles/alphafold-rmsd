import requests

uniprot = 'Q16644'

url = f'https://www.ebi.ac.uk/pdbe/graph-api/uniprot/{uniprot}'

req = requests.get(url=url)

print(req.status_code)

req_json = req.json()

print(req_json)