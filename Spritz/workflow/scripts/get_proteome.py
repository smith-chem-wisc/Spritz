# see https://www.ebi.ac.uk/training/online/sites/ebi.ac.uk.training.online/files/UniProt_programmatically_py3.pdf

import requests
import sys
import yaml

# find proteome
BASE_URL = 'https://rest.uniprot.org'
ENDPOINT = '/proteomes/search'

params = {
    'query': '*',
    'format': 'tsv',
}

with open("config/config.yaml", 'r') as stream:
   data = yaml.safe_load(stream)
organism = data["organism"].lower()

proteome_res = requests.get(BASE_URL + ENDPOINT, params=params, stream=True)
proteome_res.raise_for_status() # throw an error for bad status code

results = proteome_res.text.split('\n')[1:]
proteome = None
for r in results:
    splt = r.split('\t')
    if organism.replace('_', ' ') in splt[1].lower():
        proteome = splt[0]
        break

if proteome is None:
    print(f"Proteome for organism {organism} not found.")
    sys.exit(1)