"""Seam A: bin/write_samplesheet_from_fastq_dir.py writes one NVD row per sample, lanes as globs."""

import json
import subprocess
import sys
from pathlib import Path

import polars as pl

WRITE_SCRIPT = Path(__file__).parents[2] / "bin" / "write_samplesheet_from_fastq_dir.py"
RESOLVE_SCRIPT = Path(__file__).parents[2] / "bin" / "resolve_read_inputs.py"


def touch(directory: Path, *names: str) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for name in names:
        (directory / name).write_bytes(b"")


def write_samplesheet(tmp_path: Path, fastq_dir: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(WRITE_SCRIPT), "--fastq-dir", str(fastq_dir), "--out", str(tmp_path / "samplesheet.csv")],
        capture_output=True, text=True,
    )


def test_one_row_per_sample_with_lane_globs_that_resolve_in_lane_order(tmp_path):
    fastq_dir = tmp_path / "fastqs"
    touch(fastq_dir, *[f"{sample}_L00{lane}_R{read}_001.fastq.gz"
                       for sample in ("MO_Columbia_RWWTP_20261006_S18", "TX_Austin_20261006_S3")
                       for lane in (1, 2, 3) for read in (1, 2)])

    completed = write_samplesheet(tmp_path, fastq_dir)
    assert completed.returncode == 0, completed.stderr
    samplesheet = pl.read_csv(tmp_path / "samplesheet.csv", infer_schema=False)
    assert samplesheet.columns == ["sample_id", "srr", "platform", "fastq1", "fastq2", "fastq1_glob", "fastq2_glob"]
    assert samplesheet["sample_id"].to_list() == ["MO_Columbia_RWWTP_20261006", "TX_Austin_20261006"]
    assert samplesheet["platform"].to_list() == ["illumina", "illumina"]

    # NVD's own resolver turns each row into its lanes, in order, R1 matched with R2.
    subprocess.run([sys.executable, str(RESOLVE_SCRIPT), "--samplesheet", str(tmp_path / "samplesheet.csv"),
                    "--output-jsonl", str(tmp_path / "resolved.jsonl")], check=True)
    columbia = json.loads((tmp_path / "resolved.jsonl").read_text().splitlines()[0])
    assert [Path(path).name for path in columbia["r1"]] == [f"MO_Columbia_RWWTP_20261006_S18_L00{lane}_R1_001.fastq.gz" for lane in (1, 2, 3)]
    assert [Path(path).name for path in columbia["r2"]] == [f"MO_Columbia_RWWTP_20261006_S18_L00{lane}_R2_001.fastq.gz" for lane in (1, 2, 3)]


def test_a_sample_missing_a_mate_or_a_lane_is_refused(tmp_path):
    fastq_dir = tmp_path / "fastqs"
    touch(fastq_dir, "A_S1_L001_R1_001.fastq.gz", "A_S1_L001_R2_001.fastq.gz", "A_S1_L002_R1_001.fastq.gz")

    completed = write_samplesheet(tmp_path, fastq_dir)
    assert completed.returncode == 1
    assert "A_S1" in completed.stderr and "L002" in completed.stderr


def test_an_unfinished_rsync_temp_file_is_refused(tmp_path):
    fastq_dir = tmp_path / "fastqs"
    touch(fastq_dir, "A_S1_L001_R1_001.fastq.gz", "A_S1_L001_R2_001.fastq.gz", ".A_S1_L002_R1_001.fastq.gz.Xa91kq")

    completed = write_samplesheet(tmp_path, fastq_dir)
    assert completed.returncode == 1
    assert "transfer" in completed.stderr.lower()
