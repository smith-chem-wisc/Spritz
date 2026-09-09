# see https://www.ebi.ac.uk/training/online/sites/ebi.ac.uk.training.online/files/UniProt_programmatically_py3.pdf
#
# Resolves the UniProt proteome to download. Imported by download_uniprot.py, which reads the
# module-level `proteome` set here, so this runs at import time.

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
division = (data.get("division") or "vertebrates").lower()

# species_EnsemblBacteria.txt, downloaded by download_ensembl_bacteria_metadata, which carries the
# taxonomy id the bacterial lookup below needs. Parsed by the sibling module that owns the rest of
# this file's layout, rather than re-derived here.
BACTERIA_METADATA = "../resources/ensembl/species_EnsemblBacteria.txt"


def taxonomy_id(species):
    """The NCBI taxonomy id for an Ensembl Bacteria species directory name."""
    import ensembl_bacteria  # a sibling in scripts/, which sys.path[0] covers when run as a script

    with open(BACTERIA_METADATA, encoding="utf-8", errors="replace") as handle:
        return ensembl_bacteria.parse_taxonomy_ids(handle.read()).get(species)


def search(query):
    """The first proteome id a query returns, or None."""
    response = requests.get(BASE_URL + ENDPOINT,
                            params={'query': query, 'format': 'tsv'}, stream=True)
    response.raise_for_status()
    for line in response.text.split('\n')[1:]:
        if line.strip():
            return line.split('\t')[0]
    return None


def bacterial_proteome(species):
    """Look a bacterial proteome up by taxonomy id rather than by name.

    The name-scanning path below cannot reach these. Its query is '*' with no size parameter, so
    UniProt returns its default first page - 25 popular organisms, which is why human, mouse and
    yeast resolve - and no bacterium outside that list is ever seen. Nor does searching by name work
    reliably: "Pseudomonas aeruginosa PAO1" matches the PAO1-GFP, PAO1-VE13 and PAO1-VE2 derivatives
    before PAO1 itself, whose UniProt name is the strain-designation form
    "Pseudomonas aeruginosa (strain ATCC 15692 / ... / PAO1)".

    The taxonomy id is exact, and Ensembl Bacteria already publishes it for every species.
    """
    taxon = taxonomy_id(species)
    if taxon is None:
        print(f"No taxonomy id for {species} in {BACTERIA_METADATA}.", file=sys.stderr)
        return None

    # reference:true picks the reference proteome, of which there is at most one. PAO1 has six
    # proteomes for taxon 208964 and only UP000002438 is the reference; taking whichever came first
    # would sometimes take a redundant assembly of the same organism.
    found = search(f'organism_id:{taxon} AND reference:true')
    if found is None:
        found = search(f'organism_id:{taxon}')
        if found is not None:
            print(f"UniProt has no reference proteome for taxon {taxon}; using {found}.",
                  file=sys.stderr)
    return found


if division == "bacteria":
    proteome = bacterial_proteome(data["species"].lower())
else:
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
