"""Tests that envs/smoke-commands.tsv still matches the workflow it is meant to cover.

conda-envs.yml builds every environment in envs/ and then runs the commands this manifest lists
for it, because an environment can solve perfectly and still not provide the executable a rule
calls - see issue #237, where bin/TransDecoder.LongOrfs was a symlink to a file the package had
stopped shipping.

The manifest is maintained by hand, so the risk is that it drifts: a rule starts calling a new
tool, or a rule's `conda:` directive moves to an environment nobody listed commands for, and the
CI job goes on passing while covering less than it appears to. These checks are offline and fast,
and run in the existing workflow-script test job rather than in the slow environment job.
"""
import os
import re

import pytest

WORKFLOW_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
ENVS_DIR = os.path.join(WORKFLOW_DIR, "envs")
MANIFEST = os.path.join(ENVS_DIR, "smoke-commands.tsv")

# `conda: "../envs/x.yaml"` in rules/*.smk, `conda: "envs/x.yaml"` in the Snakefile, sometimes with
# a trailing comment. Only the basename matters here.
CONDA_DIRECTIVE = re.compile(r'^\s*conda:\s*"([^"]+)"')


def snakefiles():
    """The Snakefile and every rules/*.smk, as absolute paths."""
    paths = [os.path.join(WORKFLOW_DIR, "Snakefile")]
    rules_dir = os.path.join(WORKFLOW_DIR, "rules")
    paths += [
        os.path.join(rules_dir, name)
        for name in sorted(os.listdir(rules_dir))
        if name.endswith(".smk")
    ]
    return paths


def envs_used_by_rules():
    """Basenames of the environment files that `conda:` directives point at."""
    used = set()
    for path in snakefiles():
        with open(path) as handle:
            for line in handle:
                match = CONDA_DIRECTIVE.match(line)
                if match:
                    used.add(os.path.basename(match.group(1)))
    return used


def manifest_rows():
    """(env file, command) pairs from the manifest, ignoring comments and blank lines."""
    rows = []
    with open(MANIFEST) as handle:
        for number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            assert len(fields) == 2, f"{MANIFEST}:{number}: expected two tab-separated fields"
            rows.append((fields[0], fields[1]))
    return rows


def env_files_on_disk():
    return {name for name in os.listdir(ENVS_DIR) if name.endswith(".yaml")}


def test_manifest_names_only_environments_that_exist():
    """A typo'd environment name would silently cover nothing, since the job selects by name."""
    listed = {env for env, _ in manifest_rows()}
    missing = sorted(listed - env_files_on_disk())
    assert not missing, f"manifest names environment files that do not exist: {missing}"


def test_every_rule_environment_has_at_least_one_command():
    """The job builds an unlisted environment but checks nothing in it, so this is the real gap.

    tests.yaml is deliberately exempt: no rule uses it, it is the environment for these tests.
    """
    listed = {env for env, _ in manifest_rows()}
    uncovered = sorted(envs_used_by_rules() - listed)
    assert not uncovered, (
        "these environments are used by a rule but have no smoke commands, so conda-envs.yml "
        f"would build them and check nothing: {uncovered}"
    )


def test_commands_are_not_empty():
    for env, command in manifest_rows():
        assert command.strip(), f"{env} has an empty command"


@pytest.mark.parametrize("env", sorted(envs_used_by_rules()))
def test_rule_environments_exist_on_disk(env):
    """A `conda:` directive pointing at a missing file fails only when that rule finally runs."""
    assert env in env_files_on_disk(), f"a rule references envs/{env}, which does not exist"
