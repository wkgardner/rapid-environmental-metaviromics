"""Guards for conf/ospool.config: mistakes here only show up after a job reaches the cluster."""

import re
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OSPOOL_CONFIG = (PROJECT_ROOT / "conf" / "ospool.config").read_text()


def test_every_process_has_a_container():
    """A process without a container would look for its tools on an execute node and fail."""
    selectors = [
        selector
        for block_selector, block in re.findall(r"withName: '([^']+)' \{(.*?)\n    \}", OSPOOL_CONFIG, re.DOTALL)
        if "container =" in block
        for selector in block_selector.split("|")
    ]
    process_names = [
        re.search(r"^process (\w+)", module_file.read_text(), re.MULTILINE).group(1)
        for module_file in (PROJECT_ROOT / "modules" / "local").glob("*/main.nf")
    ]
    assert process_names
    uncovered = [
        name for name in process_names if not any(re.fullmatch(selector, name) for selector in selectors)
    ]
    assert uncovered == []


def test_auto_staging_never_copies_the_whole_project():
    """Staging projectDir itself would copy real-samples/, runs/ and ref/ onto /staging."""
    assert 'autoStageDirectories = ["${projectDir}/assets"]' in OSPOOL_CONFIG


def test_jobs_are_limited_to_staging_nodes():
    """A bare HasCHTCStaging attribute has no effect on matching, and the OS check needs parentheses."""
    assert "((OpSysMajorVer == 8) || (OpSysMajorVer == 9)) && (Target.HasCHTCStaging == true)" in OSPOOL_CONFIG
