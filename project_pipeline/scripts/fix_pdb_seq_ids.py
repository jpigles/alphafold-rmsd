import json
import requests

uniprot = 'Q16644'

url = f'https://www.ebi.ac.uk/pdbe/graph-api/uniprot/:{uniprot}'

req = requests.post(url=url)

result = json.loads(req)

print(result)