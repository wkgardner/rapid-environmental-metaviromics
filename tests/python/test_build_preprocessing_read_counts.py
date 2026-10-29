"""Seam A: bin/build_preprocessing_read_counts.py turns Deacon and BBTools logs into one table row.

Expected numbers are read off the real fixture logs by eye (tests/python/fixtures/read_qc/).
"""

import subprocess
import sys
from pathlib import Path

import polars as pl

BUILD_SCRIPT = Path(__file__).parents[2] / "bin" / "build_preprocessing_read_counts.py"
FIXTURES = Path(__file__).parent / "fixtures" / "read_qc"
COLUMNS = [
    "run_id", "sample_id", "raw_reads", "reads_removed_by_reference_filter", "reads_post_reference_filter",
    "qc_adapter_trimmed_reads", "qc_duplicate_reads_removed", "qc_low_quality_reads_removed",
    "qc_low_complexity_reads_removed", "reads_post_qc",
]


def build(tmp_path: Path, sample_id: str, *step_arguments: str) -> dict:
    completed = subprocess.run(
        [sys.executable, str(BUILD_SCRIPT), "--run-id", "rww_20261029_test", "--sample-id", sample_id,
         "--out-dir", str(tmp_path), *step_arguments],
        capture_output=True, text=True,
    )
    assert completed.returncode == 0, completed.stderr
    table_path = tmp_path / f"rww_20261029_test_{sample_id}.preprocessing_read_counts.tsv"
    table = pl.read_csv(table_path, separator="\t", infer_schema=False)
    assert table.columns == COLUMNS
    assert table.height == 1
    return table.row(0, named=True)


def test_every_step_ran(tmp_path):
    row = build(
        tmp_path, "both_callers",
        "--reference-filter-summary", str(FIXTURES / "both_callers.reference_filter.deacon.json"),
        "--adapter-trimming-log", str(FIXTURES / "both_callers.adapter_trimming.bbduk.log"),
        "--dedup-log", str(FIXTURES / "both_callers.dedup.clumpify.log"),
        "--read-filtering-log", str(FIXTURES / "both_callers.read_filtering.bbduk.log"),
    )
    assert row == {
        "run_id": "rww_20261029_test", "sample_id": "both_callers",
        "raw_reads": "4400", "reads_removed_by_reference_filter": "600", "reads_post_reference_filter": "3800",
        "qc_adapter_trimmed_reads": "246", "qc_duplicate_reads_removed": "120",
        "qc_low_quality_reads_removed": "200", "qc_low_complexity_reads_removed": "16",
        "reads_post_qc": "3584",
    }


def test_skipped_steps_read_skipped_and_reads_post_qc_comes_from_the_last_step_that_ran(tmp_path):
    row = build(
        tmp_path, "both_callers",
        "--reference-filter-summary", "skipped",
        "--adapter-trimming-log", str(FIXTURES / "both_callers.adapter_trimming.bbduk.log"),
        "--dedup-log", str(FIXTURES / "both_callers.dedup.clumpify.log"),
        "--read-filtering-log", "skipped",
    )
    # With the reference filter skipped, the raw count is what the first QC step read in.
    assert row["raw_reads"] == "3800"
    assert row["reads_removed_by_reference_filter"] == row["reads_post_reference_filter"] == "skipped"
    assert row["qc_low_quality_reads_removed"] == row["qc_low_complexity_reads_removed"] == "skipped"
    assert row["reads_post_qc"] == "3680"


def test_sample_with_no_reads_left_after_the_reference_filter_has_blank_qc_columns(tmp_path):
    row = build(tmp_path, "host_only",
                "--reference-filter-summary", str(FIXTURES / "host_only.reference_filter.deacon.json"))
    assert (row["raw_reads"], row["reads_removed_by_reference_filter"], row["reads_post_reference_filter"]) == ("600", "600", "0")
    assert row["qc_adapter_trimmed_reads"] is None
    assert row["qc_duplicate_reads_removed"] is None
    assert row["reads_post_qc"] == "0"


def test_ont_sample_without_adapter_trimming_has_a_blank_trim_column(tmp_path):
    row = build(
        tmp_path, "ont_sample",
        "--reference-filter-summary", str(FIXTURES / "both_callers.reference_filter.deacon.json"),
        "--dedup-log", str(FIXTURES / "both_callers.dedup.clumpify.log"),
        "--read-filtering-log", str(FIXTURES / "both_callers.read_filtering.bbduk.log"),
    )
    assert row["qc_adapter_trimmed_reads"] is None
    assert row["qc_duplicate_reads_removed"] == "120"
