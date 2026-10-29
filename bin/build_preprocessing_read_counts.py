#!/usr/bin/env python3
"""Build one sample's preprocessing_read_counts row from the Deacon summary and BBTools logs.

Each step argument takes one of:
    <path>     the step ran; its counts are read from this file
    skipped    the step was skipped by a --skip_<process> param; its columns read "skipped"
    (absent)   the step did not run for this sample (no reads left, or adapter trimming on
               ONT reads); its columns are left empty, which LabKey reads as null

Read counts are reads (mates), not pairs. BBDuk counts reads dropped by minlength inside
"Low quality discards", so qc_low_quality_reads_removed covers quality and length together.

Output: <out-dir>/<run_id>_<sample_id>.preprocessing_read_counts.tsv
"""

import argparse
import json
import re
from pathlib import Path

import polars as pl

SKIPPED = "skipped"
COLUMNS = [
    "run_id", "sample_id", "raw_reads", "host_reads_removed", "reads_post_depletion",
    "qc_adapter_trimmed_reads", "qc_duplicate_reads_removed", "qc_low_quality_reads_removed",
    "qc_low_complexity_reads_removed", "reads_post_qc",
]


def read_log_count(log_text: str, line_label: str) -> int | None:
    """Return the first integer after a BBTools stats label, e.g. 'KTrimmed:   246 reads' -> 246."""
    match = re.search(rf"^{re.escape(line_label)}\s+(\d+)", log_text, flags=re.MULTILINE)
    return int(match.group(1)) if match else None


def parse_bbduk_log(log_text: str) -> dict[str, int | None]:
    return {
        "reads_in": read_log_count(log_text, "Input:"),
        "reads_out": read_log_count(log_text, "Result:"),
        "adapter_trimmed_reads": read_log_count(log_text, "KTrimmed:"),
        "low_quality_reads_removed": read_log_count(log_text, "Low quality discards:"),
        "low_complexity_reads_removed": read_log_count(log_text, "Low entropy discards:"),
    }


def parse_clumpify_log(log_text: str) -> dict[str, int | None]:
    return {
        "reads_in": read_log_count(log_text, "Reads In:"),
        "reads_out": read_log_count(log_text, "Reads Out:"),
        "duplicate_reads_removed": read_log_count(log_text, "Duplicates Found:"),
    }


def load_step(argument: str | None, parse_log) -> dict | str | None:
    """A step's parsed counts, the string 'skipped', or None when the step did not run."""
    if argument is None or argument == SKIPPED:
        return argument
    return parse_log(Path(argument).read_text())


def as_cell(value: int | str | None) -> str | None:
    return None if value is None else str(value)


def build_preprocessing_read_counts_row(run_id: str, sample_id: str, host_depletion_summary: str | None,
                                        adapter_trimming_log: str | None, dedup_log: str | None,
                                        read_filtering_log: str | None) -> dict[str, str | None]:
    host_depletion = load_step(host_depletion_summary, lambda text: json.loads(text))
    adapter_trimming = load_step(adapter_trimming_log, parse_bbduk_log)
    dedup = load_step(dedup_log, parse_clumpify_log)
    read_filtering = load_step(read_filtering_log, parse_bbduk_log)
    qc_steps_in_order = [adapter_trimming, dedup, read_filtering]
    qc_steps_that_ran = [step for step in qc_steps_in_order if isinstance(step, dict)]

    def step_value(step, key):
        return step if step in (None, SKIPPED) else step[key]

    if isinstance(host_depletion, dict):
        raw_reads = host_depletion["seqs_in"]
        reads_post_depletion = host_depletion["seqs_out"]
    else:
        # Depletion was skipped: the raw count is what the first QC step read in.
        raw_reads = qc_steps_that_ran[0]["reads_in"] if qc_steps_that_ran else None
        reads_post_depletion = host_depletion

    # The reads both callers see: the output of the last QC step that ran.
    reads_post_qc = qc_steps_that_ran[-1]["reads_out"] if qc_steps_that_ran else (
        reads_post_depletion if reads_post_depletion != SKIPPED else raw_reads)

    return {
        "run_id": run_id,
        "sample_id": sample_id,
        "raw_reads": as_cell(raw_reads),
        "host_reads_removed": as_cell(step_value(host_depletion, "seqs_removed")),
        "reads_post_depletion": as_cell(reads_post_depletion),
        "qc_adapter_trimmed_reads": as_cell(step_value(adapter_trimming, "adapter_trimmed_reads")),
        "qc_duplicate_reads_removed": as_cell(step_value(dedup, "duplicate_reads_removed")),
        "qc_low_quality_reads_removed": as_cell(step_value(read_filtering, "low_quality_reads_removed")),
        "qc_low_complexity_reads_removed": as_cell(step_value(read_filtering, "low_complexity_reads_removed")),
        "reads_post_qc": as_cell(reads_post_qc),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--host-depletion-summary", help="Deacon --summary JSON, or 'skipped'")
    parser.add_argument("--adapter-trimming-log", help="BBDuk adapter-trimming stats log, or 'skipped'")
    parser.add_argument("--dedup-log", help="Clumpify stats log, or 'skipped'")
    parser.add_argument("--read-filtering-log", help="BBDuk read-filtering stats log, or 'skipped'")
    arguments = parser.parse_args()

    row = build_preprocessing_read_counts_row(
        arguments.run_id, arguments.sample_id, arguments.host_depletion_summary,
        arguments.adapter_trimming_log, arguments.dedup_log, arguments.read_filtering_log)
    arguments.out_dir.mkdir(parents=True, exist_ok=True)
    output_path = arguments.out_dir / f"{arguments.run_id}_{arguments.sample_id}.preprocessing_read_counts.tsv"
    pl.DataFrame([row], schema={column: pl.String for column in COLUMNS}).write_csv(
        output_path, separator="\t", null_value="")


if __name__ == "__main__":
    main()
