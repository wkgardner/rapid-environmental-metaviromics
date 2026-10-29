#!/usr/bin/env python3
"""Write one sample x step row of the sample_status table.

The status comes from exactly one of:
    --status <status>        given by the step itself (e.g. the EsViritu wrapper, or "skipped")
    --reads-left <count>     preprocessing: ok, or no_reads_left when no read survived
    --calls-table <path>     a call table: ok when it has calls, no_detections when header-only

Output: <out-dir>/<run_id>_<sample_id>_<step>.sample_status.tsv (one file per step, written
once when that step ends, because LabKey's file drop cannot update rows).
"""

import argparse
from pathlib import Path

import polars as pl

STEPS = ["preprocessing", "genome_kmer_based_calling", "genome_mapping_based_calling"]
STATUSES = ["ok", "no_reads_left", "skipped_no_reads", "no_detections", "failed", "skipped"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--step", required=True, choices=STEPS)
    status_source = parser.add_mutually_exclusive_group(required=True)
    status_source.add_argument("--status", choices=STATUSES)
    status_source.add_argument("--reads-left", type=int)
    status_source.add_argument("--calls-table", type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    arguments = parser.parse_args()

    if arguments.status is not None:
        status = arguments.status
    elif arguments.reads_left is not None:
        status = "ok" if arguments.reads_left > 0 else "no_reads_left"
    else:
        calls = pl.read_csv(arguments.calls_table, separator="\t", infer_schema=False)
        status = "ok" if calls.height > 0 else "no_detections"

    arguments.out_dir.mkdir(parents=True, exist_ok=True)
    pl.DataFrame([{"run_id": arguments.run_id, "sample_id": arguments.sample_id,
                   "step": arguments.step, "status": status}]).write_csv(
        arguments.out_dir / f"{arguments.run_id}_{arguments.sample_id}_{arguments.step}.sample_status.tsv",
        separator="\t")


if __name__ == "__main__":
    main()
