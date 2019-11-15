# see https://www.ebi.ac.uk/training/online/sites/ebi.ac.uk.training.online/files/UniProt_programmatically_py3.pdf

import requests
import sys

BASE = 'http://www.uniprot.org'
KB_ENDPOINT = '/uniprot/'
TOOL_ENDPOINT = '/uploadlists/'

# query = 'name:"polymerase alpha" AND proteome:UP000005640 AND reviewed:yes'
# query = 'proteome:UP000005640 AND reviewed:yes'
query = 'proteome:UP000005640'

payload = {
    'query': query,
    'format': 'xml',
    # 'columns': 'id,entry_name,reviewed,protein_names,organism,ec,keywords',
    }

result = requests.get(BASE + KB_ENDPOINT, params=payload, stream=True)
result.raise_for_status() # throw an error for bad status code
for block in result.iter_content(1024):
    sys.stdout.buffer.write(block)
