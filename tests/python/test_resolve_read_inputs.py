"""Seam A: bin/resolve_read_inputs.py turns an NVD samplesheet into one JSON record per sample."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

RESOLVE_SCRIPT = Path(__file__).parents[2] / "bin" / "resolve_read_inputs.py"
HEADER = "sample_id,srr,platform,fastq1,fastq2,fastq1_glob,fastq2_glob\n"


def touch_fastqs(directory: Path, *names: str) -> list[Path]:
    directory.mkdir(parents=True, exist_ok=True)
    paths = [directory / name for name in names]
    for path in paths:
        path.write_bytes(b"")
    return paths


def resolve(tmp_path: Path, samplesheet_rows: str) -> subprocess.CompletedProcess:
    samplesheet = tmp_path / "samplesheet.csv"
    samplesheet.write_text(HEADER + samplesheet_rows)
    return subprocess.run(
        [sys.executable, str(RESOLVE_SCRIPT), "--samplesheet", str(samplesheet),
         "--output-jsonl", str(tmp_path / "resolved_reads.jsonl")],
        capture_output=True, text=True,
    )


def resolved_records(tmp_path: Path) -> list[dict]:
    lines = (tmp_path / "resolved_reads.jsonl").read_text().splitlines()
    return [json.loads(line) for line in lines]


def test_local_paired_fastqs_resolve_to_one_paired_record(tmp_path):
    read1, read2 = touch_fastqs(tmp_path / "reads", "ww_madison_R1.fastq.gz", "ww_madison_R2.fastq.gz")
    completed = resolve(tmp_path, f"ww_madison,,illumina,{read1},{read2},,\n")

    assert completed.returncode == 0, completed.stderr
    [record] = resolved_records(tmp_path)
    assert record["sample_id"] == "ww_madison"
    assert record["source"] == "paired_files"
    assert record["read_structure"] == "paired"
    assert record["r1"] == [str(read1)] and record["r2"] == [str(read2)]


def test_sra_accession_row_resolves_to_an_sra_record(tmp_path):
    completed = resolve(tmp_path, "sra_sample,SRR33296246,illumina,,,,\n")

    assert completed.returncode == 0, completed.stderr
    [record] = resolved_records(tmp_path)
    assert record["source"] == "sra"
    assert record["srr"] == "SRR33296246"


def test_lane_split_globs_resolve_to_ordered_mate_lists(tmp_path):
    lane_files = touch_fastqs(
        tmp_path / "lanes",
        "s1_S7_L002_R1_001.fastq.gz", "s1_S7_L001_R1_001.fastq.gz",
        "s1_S7_L001_R2_001.fastq.gz", "s1_S7_L002_R2_001.fastq.gz",
    )
    lanes = lane_files[0].parent
    completed = resolve(tmp_path, f"s1,,illumina,,,{lanes}/s1_S7_L*_R1_001.fastq.gz,{lanes}/s1_S7_L*_R2_001.fastq.gz\n")

    assert completed.returncode == 0, completed.stderr
    [record] = resolved_records(tmp_path)
    assert record["source"] == "paired_globs"
    assert [Path(path).name for path in record["r1"]] == ["s1_S7_L001_R1_001.fastq.gz", "s1_S7_L002_R1_001.fastq.gz"]
    assert [Path(path).name for path in record["r2"]] == ["s1_S7_L001_R2_001.fastq.gz", "s1_S7_L002_R2_001.fastq.gz"]


def test_ont_barcode_folder_glob_resolves_to_single_reads(tmp_path):
    barcode_files = touch_fastqs(tmp_path / "barcode01", "chunk_0.fastq.gz", "chunk_1.fastq.gz")
    completed = resolve(tmp_path, f"barcode01,,nanopore,,,{barcode_files[0].parent}/*.fastq.gz,\n")

    assert completed.returncode == 0, completed.stderr
    [record] = resolved_records(tmp_path)
    assert record["platform"] == "ont"
    assert record["source"] == "single_glob"
    assert record["read_structure"] == "single"
    assert len(record["reads"]) == 2


@pytest.mark.parametrize("rows", [
    "dup,SRR1,illumina,,,,\ndup,SRR2,illumina,,,,\n",
    "bad id,SRR1,illumina,,,,\n",
])
def test_duplicate_or_unsafe_sample_ids_are_rejected(tmp_path, rows):
    completed = resolve(tmp_path, rows)

    assert completed.returncode == 1
    assert "sample" in completed.stderr.lower()
