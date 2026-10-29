#!/usr/bin/env python3
"""Write an NVD-format samplesheet from a folder of Illumina FASTQs, one row per sample.

Files must be named <sample>_S<n>_L<lane>_R<1|2>_001.fastq.gz (bcl2fastq / BCL Convert).
Each sample gets one row whose fastq1_glob / fastq2_glob match all of its lanes, so NVD's
resolver (bin/resolve_read_inputs.py) pairs the lanes in order. sample_id drops the _S<n>
sample number, e.g. MO_Columbia_RWWTP_20261006_S18 -> MO_Columbia_RWWTP_20261006.

The folder is refused when a transfer looks unfinished (rsync's hidden .<name>.XXXXXX temp
files) or when any sample lacks a mate for a lane, or a lane present for other files.
"""

import argparse
import re
import sys
from pathlib import Path

import polars as pl

FASTQ_NAME = re.compile(r"^(?P<prefix>(?P<sample>.+)_S\d+)_L(?P<lane>\d{3})_R(?P<read>[12])_001\.fastq\.gz$")
NVD_COLUMNS = ["sample_id", "srr", "platform", "fastq1", "fastq2", "fastq1_glob", "fastq2_glob"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fastq-dir", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--platform", default="illumina")
    arguments = parser.parse_args()

    fastq_dir = arguments.fastq_dir.resolve()
    unfinished = sorted(path.name for path in fastq_dir.iterdir() if path.name.startswith(".") and ".fastq.gz" in path.name)
    if unfinished:
        sys.stderr.write(f"Transfer not finished, rsync temp files present: {', '.join(unfinished)}\n")
        raise SystemExit(1)

    lanes_by_prefix: dict[str, dict[str, set[str]]] = {}
    for path in sorted(fastq_dir.iterdir()):
        match = FASTQ_NAME.match(path.name)
        if match:
            reads = lanes_by_prefix.setdefault(match["prefix"], {"1": set(), "2": set()})
            reads[match["read"]].add(match["lane"])

    problems = []
    for prefix, reads in sorted(lanes_by_prefix.items()):
        for lane in sorted(reads["1"] ^ reads["2"]):
            problems.append(f"{prefix}: lane L{lane} has only one mate")
    all_lanes = set().union(*(reads["1"] | reads["2"] for reads in lanes_by_prefix.values())) if lanes_by_prefix else set()
    for prefix, reads in sorted(lanes_by_prefix.items()):
        for lane in sorted(all_lanes - (reads["1"] & reads["2"]) - (reads["1"] ^ reads["2"])):
            problems.append(f"{prefix}: lane L{lane} is missing")
    if problems:
        sys.stderr.write("Incomplete samples:\n" + "\n".join(problems) + "\n")
        raise SystemExit(1)

    rows = [{
        "sample_id": FASTQ_NAME.match(f"{prefix}_L001_R1_001.fastq.gz")["sample"],
        "srr": "", "platform": arguments.platform, "fastq1": "", "fastq2": "",
        "fastq1_glob": f"{fastq_dir}/{prefix}_L[0-9][0-9][0-9]_R1_001.fastq.gz",
        "fastq2_glob": f"{fastq_dir}/{prefix}_L[0-9][0-9][0-9]_R2_001.fastq.gz",
    } for prefix in sorted(lanes_by_prefix)]
    pl.DataFrame(rows, schema={column: pl.String for column in NVD_COLUMNS}).write_csv(arguments.out)
    print(f"{len(rows)} samples, lanes {', '.join(f'L{lane}' for lane in sorted(all_lanes))} -> {arguments.out}")


if __name__ == "__main__":
    main()
