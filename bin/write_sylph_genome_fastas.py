#!/usr/bin/env python3
"""Write one FASTA per VPD Assembly for the sylph database, cut at N runs.

sylph encodes N as A (sylph src/types.rs), so masked or ambiguous bases would add false k-mers.
Every record is cut at every run of N, and fragments shorter than --min-fragment-length (the
k-mer length) are dropped; fragments are named <accession>__f<index>. All of an Assembly's
records (segments) go in its one file, so a segmented virus is one sylph genome.

Outputs: <out-dir>/<genome_file> per Assembly (names from the VPD lookup), and --genome-list,
one path per line, for `sylph sketch --gl`.
"""

import argparse
import re
from pathlib import Path

import polars as pl
from Bio import SeqIO


def cut_at_n_runs(sequence: str, min_fragment_length: int) -> list[str]:
    return [fragment for fragment in re.split(r"[Nn]+", sequence) if len(fragment) >= min_fragment_length]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vpd-fasta", required=True, type=Path)
    parser.add_argument("--vpd-metadata", required=True, type=Path)
    parser.add_argument("--vpd-lookup", required=True, type=Path, help="from build_vpd_lookup.py")
    parser.add_argument("--min-fragment-length", required=True, type=int)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--genome-list", required=True, type=Path)
    arguments = parser.parse_args()

    metadata = pl.read_csv(arguments.vpd_metadata, separator="\t", infer_schema=False, quote_char=None)
    lookup = pl.read_csv(arguments.vpd_lookup, separator="\t", infer_schema=False, quote_char=None)
    genome_file_by_accession = dict(
        metadata.join(lookup, left_on="Assembly", right_on="vpd_assembly").select("Accession", "genome_file").iter_rows())

    arguments.out_dir.mkdir(parents=True, exist_ok=True)
    fragments_by_genome_file: dict[str, list[str]] = {}
    for record in SeqIO.parse(arguments.vpd_fasta, "fasta"):
        genome_file = genome_file_by_accession[record.id]
        fragments = cut_at_n_runs(str(record.seq), arguments.min_fragment_length)
        fragments_by_genome_file.setdefault(genome_file, []).extend(
            f">{record.id}__f{index}\n{fragment}\n" for index, fragment in enumerate(fragments))

    genome_paths = []
    for genome_file, fasta_records in sorted(fragments_by_genome_file.items()):
        genome_path = arguments.out_dir / genome_file
        genome_path.write_text("".join(fasta_records))
        genome_paths.append(str(genome_path))
    arguments.genome_list.write_text("\n".join(genome_paths) + "\n")


if __name__ == "__main__":
    main()
