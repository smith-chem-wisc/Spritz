"""The run's configuration, as the workflow sees it.

A script cannot reach the config by relative path. Snakemake runs from `workflow/`, where
`config/config.yaml` is the *packaged default* that ships beside the assembly and is baked into the
container image; the run's own config is written into the analysis directory and handed over with
`--configfile`. Opening the relative path therefore reads the defaults - silently, because they parse
perfectly well - so a run asking for bacteria would be told it was analysing `Homo_sapiens`.

Every script that needs a config value takes the real path from its rule, through SPRITZ_CONFIG.

This is also the single home for the adjustments the workflow makes before reading a value, so the
snakefiles and the scripts cannot drift apart about what the config means. `common.smk` calls
`normalise` on snakemake's own config dict for exactly that reason.
"""

import os

import yaml

DEFAULT_PATH = "config/config.yaml"
ENV_VAR = "SPRITZ_CONFIG"


def choose_configfile(configfiles, packaged_default):
    """Pick this run's config out of snakemake's `workflow.configfiles`.

    Not `configfiles[-1]`. Snakemake seeds that list with the paths given to `--configfile` and
    *then* appends whatever the `configfile:` directive names, and Spritz's Snakefile opens with
    `configfile: "config/config.yaml"`. So the last entry is the packaged default, and picking it
    reintroduced the very bug this module exists to prevent: every script saw Homo_sapiens, a
    bacterial run downloaded the human proteome under the bacterium's filename, and prose.txt
    described a quant run that never happened.

    Select by identity rather than position - anything that is not the packaged default was
    supplied for this run - so the answer does not depend on how snakemake orders the list.
    """
    default = os.path.realpath(str(packaged_default))
    supplied = [str(p) for p in configfiles if os.path.realpath(str(p)) != default]
    return supplied[-1] if supplied else str(packaged_default)


def normalise(config):
    """Apply the adjustments the workflow relies on, in place, and return the config."""
    # `--config vcf=x.vcf` on the command line and a hand-written `vcf: "x.vcf"` both give a plain
    # string. len() on a string is its length, so check() would pass and then iterating it would
    # yield one entry per character.
    if isinstance(config.get("vcf"), str):
        config["vcf"] = [config["vcf"]] if config["vcf"] else []

    # `reference_free: "False"` in a hand-written config, or --config reference_free=False on a
    # shell that quotes it, gives the STRING "False" - and every non-empty string is truthy, so the
    # option would read as on. Same trap as `vcf` above. Only the words that mean "off" are mapped;
    # anything else keeps normal truthiness, so a typo does not silently disable the feature.
    if isinstance(config.get("reference_free"), str):
        config["reference_free"] = \
            config["reference_free"].strip().lower() not in ("", "false", "no", "0", "off", "none")

    # Ensembl Bacteria publishes no variation at all, so the downloaded route can never apply to a
    # bacterial reference; bootstrap is simply the route bacteria take.
    if (config.get("division") or "vertebrates").lower() == "bacteria" \
            and (config.get("known_sites") or "ensembl") == "ensembl":
        config["known_sites"] = "bootstrap"

    return config


def load(path=None):
    """The run's config: the path given, else SPRITZ_CONFIG, else the packaged default."""
    with open(path or os.environ.get(ENV_VAR) or DEFAULT_PATH, encoding="utf-8") as handle:
        return normalise(yaml.safe_load(handle) or {})


def check(config, field):
    """Whether a field is present and non-empty. The same test check() makes in common.smk."""
    return field in config and config[field] is not None and len(config[field]) > 0
