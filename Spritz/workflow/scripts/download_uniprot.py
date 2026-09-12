"""Downloads the UniProt proteome that get_proteome resolved, in xml or fasta.

Two things here are not obvious, and both were silent.

`rest.uniprot.org/uniprotkb/search` PAGINATES, at 25 entries. The old www.uniprot.org endpoint this
was written against returned every hit in one response, so when UniProt migrated the API the same
call quietly started returning the first page: 25 proteins for every species, human included, which
is what the modification transfer then had to work with. `/uniprotkb/stream` is the endpoint meant
for downloads and returns the whole set, so it is used instead of following `Link: rel="next"`.

And the proteome id has to go in the `proteome:` FIELD. As a bare term it is a free-text search,
which matched 25 unrelated entries for UP000002438 and nothing at all for UP000006822.
"""

import sys

BASE_URL = "https://rest.uniprot.org"
# /uniprotkb/stream, not /uniprot/search: see the module docstring. The trailing endpoint name
# changed with the API migration too - /uniprot/ still answers, which is why this went unnoticed.
ENDPOINT = "/uniprotkb/stream"
BLOCK = 1 << 16


def queries(proteome, organism_id=None):
    """The queries to try, in order, for a resolved proteome id.

    A proteome record can exist and still have no UniProtKB entries behind it. Mesomycoplasma
    hyopneumoniae 232 (UP000006822) advertises 691 proteins and returns zero for
    `proteome:UP000006822`; its 129 Swiss-Prot entries are reachable only by organism. Falling back
    by taxonomy id is the difference between a small proteome and an empty file.
    """
    tries = [f"proteome:{proteome}"]
    if organism_id:
        tries.append(f"organism_id:{organism_id}")
    return tries


def fetch(session, query, fmt):
    """Stream one query. Returns (first_block, remaining_blocks); a falsy first block means no hits.

    Deliberately does not buffer the response: a human proteome in xml is hundreds of megabytes, so
    emptiness is judged from the first block alone rather than by reading it all in to count.
    """
    response = session.get(BASE_URL + ENDPOINT,
                           params={"query": query, "format": fmt, "includeIsoform": "true"},
                           stream=True)
    response.raise_for_status()
    blocks = response.iter_content(BLOCK)
    return next(blocks, b""), blocks


def main(argv):
    import requests

    import get_proteome  # resolves the proteome at import time, and exits if it cannot

    fmt = argv[1]
    session = requests.Session()
    attempted = queries(get_proteome.proteome, getattr(get_proteome, "organism_id", None))
    for i, query in enumerate(attempted):
        first, rest = fetch(session, query, fmt)
        if not first:
            print(f"UniProt returned nothing for '{query}'.", file=sys.stderr)
            continue
        if i:
            print(f"Using '{query}' instead.", file=sys.stderr)
        sys.stdout.buffer.write(first)
        for block in rest:
            sys.stdout.buffer.write(block)
        return 0

    print(f"UniProt returned no entries for any of: {', '.join(attempted)}.", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
