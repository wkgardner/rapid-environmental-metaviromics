#!/usr/bin/env python3
"""Build the VPD lookup: one row per VPD Assembly, shared by both callers' tables.

Columns: vpd_assembly, genome_file, vpd_species, vpd_subspecies, vpd_genus, vpd_family, vpd_taxids.

- Rank prefixes are stripped ("s__Tobamovirus capsici" -> "Tobamovirus capsici").
- vpd_subspecies and vpd_taxids are ';'-joined, sorted unique values over the Assembly's
  records (e.g. rotavirus segment genotypes). vpd_subspecies is VPD's label for the reference,
  not a strain call.
- genome_file is the sylph genome FASTA name: the Assembly with every character outside
  [A-Za-z0-9_.-] replaced by "_", plus ".fna" (123 VPD v3.2.4 Assemblies contain "/").
  The table keeps the mapping explicit, so sylph's Genome_file is joined, never un-sanitised.
  Two Assemblies that sanitise to one name are an error.
"""

import argparse
import re
import sys
from pathlib import Path

import polars as pl


def sanitise_genome_file_name(assembly: str) -> str:
    """'PSVA/S2-2__Sadwavirus_alphananas' -> 'PSVA_S2-2__Sadwavirus_alphananas.fna'."""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", assembly) + ".fna"


def strip_rank_prefix(column_name: str) -> pl.Expr:
    return pl.col(column_name).str.replace(r"^[a-z]__", "")


def build_vpd_lookup(vpd_metadata: pl.DataFrame) -> pl.DataFrame:
    return (
        vpd_metadata.group_by("Assembly", maintain_order=True)
        .agg(
            strip_rank_prefix("species").first().alias("vpd_species"),
            strip_rank_prefix("subspecies").unique().sort().str.join(";").alias("vpd_subspecies"),
            strip_rank_prefix("genus").first().alias("vpd_genus"),
            strip_rank_prefix("family").first().alias("vpd_family"),
            pl.col("TaxID").unique().sort().str.join(";").alias("vpd_taxids"),
        )
        .with_columns(
            pl.col("Assembly").map_elements(sanitise_genome_file_name, return_dtype=pl.String).alias("genome_file"))
        .select(
            pl.col("Assembly").alias("vpd_assembly"), "genome_file", "vpd_species", "vpd_subspecies",
            "vpd_genus", "vpd_family", "vpd_taxids",
        )
        .sort("vpd_assembly")
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vpd-metadata", required=True, type=Path, help="VPD all_metadata.tsv")
    parser.add_argument("--out", required=True, type=Path, help="lookup TSV to write")
    arguments = parser.parse_args()

    vpd_metadata = pl.read_csv(arguments.vpd_metadata, separator="\t", infer_schema=False, quote_char=None)
    lookup = build_vpd_lookup(vpd_metadata)
    clashing_files = lookup.group_by("genome_file").agg(pl.col("vpd_assembly")).filter(pl.col("vpd_assembly").list.len() > 1)
    if clashing_files.height:
        for row in clashing_files.iter_rows(named=True):
            sys.stderr.write(f"Assemblies {row['vpd_assembly']} all sanitise to {row['genome_file']}\n")
        raise SystemExit(1)
    lookup.write_csv(arguments.out, separator="\t", quote_style="never")


if __name__ == "__main__":
    main()
