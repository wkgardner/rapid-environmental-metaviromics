"""Seam A: bin/write_run_metadata.py records what ran, with which tool versions and settings."""

import json
import subprocess
import sys
from pathlib import Path

import yaml

WRITE_SCRIPT = Path(__file__).parents[2] / "bin" / "write_run_metadata.py"
VERSIONS = Path(__file__).parent / "fixtures" / "versions"

EFFECTIVE_PARAMS = {
    "input": "/data/samplesheet.csv",
    "reference_filter_index": "/refs/panhuman1_silva.k31w15.idx",
    "reference_filter_index_description": "panhuman-1 + SILVA SSU/LSU rRNA",
    "reference_filter_mode": "deplete",
    "sylph_database": "/refs/vpd_3.2.4.k31c50.syldb",
    "esviritu_database": "/refs/esviritu_db_v3.2.4",
    "vpd_version": "3.2.4", "vpd_mask_version": "vector_only_v1", "taxdump_version": None,
    "trim_adapters": True, "dedup": False, "skip_read_filtering": True, "sylph_minimum_ani": 90,
}


def write_metadata(tmp_path: Path) -> dict:
    params_path = tmp_path / "params.json"
    params_path.write_text(json.dumps(EFFECTIVE_PARAMS))
    completed = subprocess.run(
        [sys.executable, str(WRITE_SCRIPT), "--run-id", "rww_20261029_happy_turing",
         "--pipeline-version", "1.0.0dev", "--start-time", "2026-10-29T14:00:00Z",
         "--params-json", str(params_path), "--versions", *map(str, sorted(VERSIONS.glob("*.yml"))),
         "--out-dir", str(tmp_path)],
        capture_output=True, text=True,
    )
    assert completed.returncode == 0, completed.stderr
    return yaml.safe_load((tmp_path / "rww_20261029_happy_turing.run_metadata.yaml").read_text())


def test_run_identity_and_references_are_recorded(tmp_path):
    metadata = write_metadata(tmp_path)

    assert metadata["run_id"] == "rww_20261029_happy_turing"
    assert metadata["pipeline_version"] == "1.0.0dev"
    assert metadata["start_time"] == "2026-10-29T14:00:00Z"
    assert metadata["samplesheet"] == "/data/samplesheet.csv"
    assert metadata["references"] == {"vpd_version": "3.2.4", "vpd_mask_version": "vector_only_v1", "taxdump_version": None}


def test_each_step_records_its_tool_version_and_what_it_ran_against(tmp_path):
    steps = write_metadata(tmp_path)["steps"]

    assert steps["reference_filter"] == {
        "tool": "deacon", "version": "0.18.0", "mode": "deplete",
        "index": "/refs/panhuman1_silva.k31w15.idx", "index_contents": "panhuman-1 + SILVA SSU/LSU rRNA",
    }
    assert steps["adapter_trimming"] == {"tool": "bbmap", "version": "39.91"}
    assert steps["genome_kmer_based_calling"] == {"tool": "sylph", "version": "1.0.0", "database": "/refs/vpd_3.2.4.k31c50.syldb"}
    assert steps["genome_mapping_based_calling"] == {"tool": "esviritu", "version": "1.3.5", "database": "/refs/esviritu_db_v3.2.4"}


def test_switched_off_steps_are_marked_skipped_and_every_effective_param_is_kept(tmp_path):
    metadata = write_metadata(tmp_path)

    # dedup and adapter trimming use NVD's on/off params; every other step has skip_<step>.
    assert metadata["steps"]["dedup"] == {"skipped": True}
    assert metadata["steps"]["read_filtering"] == {"skipped": True}
    assert metadata["params"] == EFFECTIVE_PARAMS
