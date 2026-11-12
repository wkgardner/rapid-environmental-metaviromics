"""Seam A: bin/write_sample_status.py writes one sample x step row of sample_status."""

import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

WRITE_SCRIPT = Path(__file__).parents[2] / "bin" / "write_sample_status.py"
FIXTURES = Path(__file__).parent / "fixtures"


def write_status(tmp_path: Path, step: str, *status_arguments: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(WRITE_SCRIPT), "--run-id", "rww_20261029_test", "--sample-id", "s1",
         "--step", step, *status_arguments, "--out-dir", str(tmp_path)],
        capture_output=True, text=True,
    )


def status_row(tmp_path: Path, step: str) -> dict:
    table = pl.read_csv(tmp_path / f"rww_20261029_test_s1_{step}.sample_status.tsv", separator="\t", infer_schema=False)
    assert table.columns == ["run_id", "sample_id", "step", "status"]
    return table.row(0, named=True)


@pytest.mark.parametrize("reads_post_qc, expected", [("3584", "ok"), ("0", "no_reads_left")])
def test_preprocessing_status_comes_from_the_reads_left(tmp_path, reads_post_qc, expected):
    completed = write_status(tmp_path, "preprocessing", "--reads-left", reads_post_qc)
    assert completed.returncode == 0, completed.stderr
    assert status_row(tmp_path, "preprocessing") == {
        "run_id": "rww_20261029_test", "sample_id": "s1", "step": "preprocessing", "status": expected}


# The pipeline passes sylph's own profile here, not a REMV table: the step only asks whether
# the caller reported anything, and a header-only profile means it did not.
@pytest.mark.parametrize("calls_table, expected", [
    ("both_callers.sylph_profile.tsv", "ok"),
    ("mapping_only.sylph_profile.tsv", "no_detections"),
])
def test_kmer_calling_status_comes_from_whether_the_table_has_calls(tmp_path, calls_table, expected):
    completed = write_status(tmp_path, "genome_kmer_based_calling", "--calls-table", str(FIXTURES / "calls" / calls_table))
    assert completed.returncode == 0, completed.stderr
    assert status_row(tmp_path, "genome_kmer_based_calling")["status"] == expected


@pytest.mark.parametrize("status", ["ok", "no_detections", "failed", "skipped", "skipped_no_reads"])
def test_an_explicit_status_is_written_as_given(tmp_path, status):
    completed = write_status(tmp_path, "genome_mapping_based_calling", "--status", status)
    assert completed.returncode == 0, completed.stderr
    assert status_row(tmp_path, "genome_mapping_based_calling")["status"] == status


def test_unknown_steps_and_statuses_are_rejected(tmp_path):
    assert write_status(tmp_path, "species_confirmation", "--status", "ok").returncode != 0
    assert write_status(tmp_path, "preprocessing", "--status", "questionable").returncode != 0
