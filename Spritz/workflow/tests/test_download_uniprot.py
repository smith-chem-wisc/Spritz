"""Tests for download_uniprot.py's query construction.

The bug these exist for was silent and affected every species, human included. UniProt's
`/uniprotkb/search` paginates at 25; the old www.uniprot.org endpoint the script was written
against did not, so after the API migration every run's UniProt proteome was its first 25 entries.
Passing the proteome id as a bare term rather than in the `proteome:` field compounded it - a
free-text match, which returned 25 unrelated entries for one proteome and none at all for another.

Verified against the live API when written: proteome:UP000002438 gives 5564 sequences from
/uniprotkb/stream against 25 from the old call, and UP000006822 gives 0 either way but 129 by
organism id.
"""
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts")))

# noqa: E402 below - the sys.path insert above has to run first
from download_uniprot import BASE_URL, ENDPOINT, queries  # noqa: E402


def test_proteome_id_goes_in_the_proteome_field():
    """A bare id is a free-text search, which is how UP000002438 matched 25 unrelated entries."""
    assert queries("UP000002438") == ["proteome:UP000002438"]


def test_organism_id_is_offered_as_a_fallback():
    """UP000006822 advertises 691 proteins and returns none; its 129 are reachable by organism."""
    assert queries("UP000006822", 295358) == [
        "proteome:UP000006822",
        "organism_id:295358",
    ]


def test_proteome_is_tried_before_organism():
    """Organism is broader - it can pull in entries from other assemblies of the same taxon."""
    assert queries("UP1", 9606)[0].startswith("proteome:")


def test_no_fallback_without_an_organism_id():
    """The vertebrate path resolves by name and never looks a taxonomy id up."""
    for absent in (None, 0, ""):
        assert queries("UP1", absent) == ["proteome:UP1"]


def test_uses_the_streaming_endpoint():
    """/uniprotkb/search would paginate at 25 however the query is written."""
    assert ENDPOINT == "/uniprotkb/stream"
    assert BASE_URL + ENDPOINT == "https://rest.uniprot.org/uniprotkb/stream"
