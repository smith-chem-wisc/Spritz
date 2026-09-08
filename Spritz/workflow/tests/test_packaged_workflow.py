"""Tests that SpritzCMD.csproj still ships every workflow script the pipeline needs.

The csproj enumerates workflow files one `<None Include>` at a time. They live outside the project
directory, so the SDK's default glob does not pick them up and a new file is simply absent from the
build output - and the Dockerfile copies that output, so the image ships a workflow missing a file.

That failure is invisible to the linting and dry-run jobs, which run against the repository tree
rather than the built artifact. It surfaces only when the container runs, and for an imported module
it surfaces as a ModuleNotFoundError in common.smk before any rule is evaluated, so nothing runs at
all. Cheaper to assert here.
"""
import os
import re

import pytest

WORKFLOW_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SCRIPTS_DIR = os.path.join(WORKFLOW_DIR, "scripts")
CSPROJ = os.path.abspath(
    os.path.join(WORKFLOW_DIR, "..", "SpritzCMD", "SpritzCMD.csproj")
)

# <None Include="..\workflow\scripts\foo.py" Link="workflow\scripts\foo.py">, backslashed because
# MSBuild wrote it on Windows.
INCLUDED_SCRIPT = re.compile(r'Include="\.\.\\workflow\\scripts\\([^"]+)"')


def included_scripts():
    with open(CSPROJ) as handle:
        return set(INCLUDED_SCRIPT.findall(handle.read()))


def scripts_in_tree():
    return {name for name in os.listdir(SCRIPTS_DIR) if name.endswith(".py")}


def test_the_csproj_exists_where_this_test_expects():
    """A moved csproj would make every other assertion here pass vacuously."""
    assert os.path.isfile(CSPROJ)
    assert included_scripts(), "no workflow scripts matched; the Include pattern has changed"


def test_every_workflow_script_is_shipped():
    missing = sorted(scripts_in_tree() - included_scripts())
    assert not missing, (
        "these scripts exist in workflow/scripts but SpritzCMD.csproj does not copy them, so they "
        f"will be absent from the built artifact and from the container image: {missing}"
    )


def test_the_csproj_ships_nothing_that_no_longer_exists():
    """A stale entry fails the build with MSB3030 rather than at runtime, but still fails it."""
    stale = sorted(included_scripts() - scripts_in_tree())
    assert not stale, f"SpritzCMD.csproj copies scripts that are not in the tree: {stale}"


@pytest.mark.parametrize("module", ["snpeff_config.py"])
def test_modules_imported_by_the_rules_are_shipped(module):
    """Modules imported at DAG-build time are the worst case: nothing runs without them.

    A plain script is only needed by the rule that shells out to it, so its absence breaks one rule
    late in a run. common.smk imports this one, so its absence breaks the workflow before the first
    job.
    """
    assert module in scripts_in_tree()
    assert module in included_scripts()
