"""Shared fixtures for the workflow script tests.

The scripts in workflow/scripts are not importable modules: each one does its work at import
time, reading config/config.yaml and ../resources/... by relative path and communicating over
stdin, stdout and sys.argv. So they are exercised here the way the workflow actually invokes
them - as subprocesses, from a directory laid out the way they expect - rather than by
importing them. That tests the real contract, and needs no changes to the scripts themselves.
"""
import os
import subprocess
import sys
import textwrap

import pytest

SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))


@pytest.fixture
def workflow_dir(tmp_path):
    """A directory shaped like workflow/, with a sibling resources/ tree.

    Several scripts read "config/config.yaml" and "../resources/..." relative to the working
    directory, so the fixture reproduces that shape instead of patching the scripts.
    """
    resources = tmp_path / "resources"
    (resources / "ChromosomeMappings").mkdir(parents=True)
    (resources / "ensembl").mkdir(parents=True)

    workflow = tmp_path / "workflow"
    (workflow / "config").mkdir(parents=True)
    (workflow / "config" / "config.yaml").write_text(
        textwrap.dedent(
            """\
            species: "Homo_sapiens"
            genome: "GRCh38"
            release: "100"
            organism: "human"
            analyses: [quant]
            analysis_directory: [../results]
            """
        )
    )
    return workflow


def run_script(name, cwd, stdin="", args=()):
    """Run a workflow script the way a rule would, returning the CompletedProcess."""
    return subprocess.run(
        [sys.executable, os.path.join(SCRIPTS_DIR, name), *args],
        cwd=str(cwd),
        input=stdin,
        capture_output=True,
        text=True,
    )
