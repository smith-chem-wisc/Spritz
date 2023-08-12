import get_proteome
import requests
import sys

format = sys.argv[1]

proteome = get_proteome.proteome
BASE_URL = 'https://rest.uniprot.org'
ENDPOINT = '/uniprot/search'

payload = {
    'query': proteome,
    'format': format,
    'includeIsoforms': 'yes', # include isoforms in fasta
    }

result = requests.get(BASE_URL + ENDPOINT, params=payload, stream=True)
result.raise_for_status() # throw an error for bad status code
for block in result.iter_content(1024):
    sys.stdout.buffer.write(block)
