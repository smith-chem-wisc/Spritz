"""Asserts every rule that runs SpritzModifications on the UniProt xml is serialized.

mzLib's ProteinDbLoader.LoadProteinXML decompresses a .gz input to a FIXED path - temp.xml beside
the input - reads it, then deletes it. Four rules in proteogenomics.smk run that assembly on the
single shared UNIPROTXML, none of them ordered against the others, so snakemake schedules them
together and they collide on that one temp.xml:

    IOException: The process cannot access the file '/app/spritz/resources/uniprot/temp.xml'
    because it is being used by another process.  at System.IO.File.Create(String path)

It only shows up at real proteome size. While download_uniprot.py was silently fetching 25 entries
each invocation was inside the decompress-and-read window for milliseconds and they never
overlapped; at 5564 entries the window is seconds and the run dies at the last step.

The `resources: uniprot_temp=1` token serializes them. This test exists because the failure mode is
reintroduced by *adding a rule* - a fifth consumer of UNIPROTXML with no token silently races
again - which no amount of testing the four current rules would catch.
"""
import os
import re

RULES_SMK = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "rules", "proteogenomics.smk")
)

# Indented too: build_transfer_mods lives inside `if not PREBUILT_SPRITZ_MODS:`.
RULE_HEADER = re.compile(r"^[ \t]*rule[ \t]+(\w+)[ \t]*:", re.MULTILINE)

# The assembly is always invoked through the input it declares, never by a bare path.
INVOKES_ASSEMBLY = "{input.transfermods}"
TAKES_UNIPROT_XML = "unixml=UNIPROTXML"


def rule_blocks():
    """{rule name: its text}, split on the rule headers."""
    text = open(RULES_SMK, encoding="utf-8").read()
    starts = [(m.group(1), m.start()) for m in RULE_HEADER.finditer(text)]
    bounds = [s for _, s in starts[1:]] + [len(text)]
    return {name: text[start:end] for (name, start), end in zip(starts, bounds)}


def rules_loading_the_uniprot_xml():
    return {
        name: body
        for name, body in rule_blocks().items()
        if INVOKES_ASSEMBLY in body and TAKES_UNIPROT_XML in body
    }


def test_the_rules_file_parses_into_rules():
    """A moved or renamed file would make every assertion below pass vacuously."""
    assert os.path.isfile(RULES_SMK)
    blocks = rule_blocks()
    assert len(blocks) >= 6, f"only parsed {len(blocks)} rules; the header pattern has changed"
    assert "transfer_modifications_variant" in blocks


def test_the_uniprot_xml_consumers_are_found_by_content_not_by_name():
    """The markers still match, so the list below is measured rather than hardcoded."""
    found = rules_loading_the_uniprot_xml()
    assert found, (
        f"no rule contains both {INVOKES_ASSEMBLY!r} and {TAKES_UNIPROT_XML!r}; "
        "the markers have changed and this whole module is now vacuous"
    )
    assert set(found) == {
        "transfer_modifications_variant",
        "transfer_modifications_isoformvariant",
        "reference_protein_xml",
        "custom_protein_xml",
    }, (
        "the set of rules running SpritzModifications on the UniProt xml has changed: "
        f"{sorted(found)}. A new one needs `resources: uniprot_temp=1` too."
    )


def test_every_uniprot_xml_consumer_declares_the_serialization_token():
    """Without the token snakemake runs them concurrently and they collide on temp.xml."""
    missing = [
        name
        for name, body in rules_loading_the_uniprot_xml().items()
        if "uniprot_temp" not in body
    ]
    assert not missing, (
        f"rules run SpritzModifications on the shared UniProt xml without "
        f"`resources: uniprot_temp=1`: {sorted(missing)}"
    )


def test_the_token_is_only_on_rules_that_need_it():
    """A token on an unrelated rule serializes it against these for no reason."""
    consumers = set(rules_loading_the_uniprot_xml())
    stray = [
        name
        for name, body in rule_blocks().items()
        if "uniprot_temp" in body and name not in consumers
    ]
    assert not stray, f"uniprot_temp on rules that do not load the UniProt xml: {sorted(stray)}"
