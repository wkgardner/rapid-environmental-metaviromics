#!/usr/bin/env python3
"""Build one sample's genome_kmer_based_calls table from its sylph profile.

Every sylph profile row is a genome k-mer-based call: one row per sample x VPD Assembly.
Rows join to the VPD lookup on the basename of sylph's Genome_file (never Contig_name, which
holds only the first segment). A Genome_file missing from the lookup means the sylph database
and the lookup come from different VPD builds, so the script stops.

Output: <out-dir>/<run_id>_<sample_id>.genome_kmer_based_calls.tsv (header only when sylph
called nothing).
"""

import argparse
import sys
from pathlib import Path

import polars as pl

LOOKUP_COLUMNS = ["vpd_assembly", "vpd_species", "vpd_subspecies", "vpd_genus", "vpd_family", "vpd_taxids"]


def build_genome_kmer_based_calls(run_id: str, sample_id: str, sylph_profile: pl.DataFrame,
                                  vpd_lookup: pl.DataFrame) -> pl.DataFrame:
    containment_parts = pl.col("Containment_ind").str.split("/")
    return (
        sylph_profile
        .with_columns(pl.col("Genome_file").str.split("/").list.last().alias("genome_file"))
        .join(vpd_lookup, on="genome_file", how="left")
        .select(
            pl.lit(run_id).alias("run_id"),
            pl.lit(sample_id).alias("sample_id"),
            *LOOKUP_COLUMNS,
            pl.col("Taxonomic_abundance").cast(pl.Float64).alias("taxonomic_abundance"),
            pl.col("Sequence_abundance").cast(pl.Float64).alias("sequence_abundance"),
            pl.col("Adjusted_ANI").cast(pl.Float64).alias("adjusted_ani"),
            pl.col("Naive_ANI").cast(pl.Float64).alias("naive_ani"),
            pl.col("Eff_cov").cast(pl.Float64).alias("effective_coverage"),
            containment_parts.list.get(0).cast(pl.Int64).alias("containment_matching_kmers"),
            containment_parts.list.get(1).cast(pl.Int64).alias("containment_total_kmers"),
            pl.col("kmers_reassigned").cast(pl.Int64),
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--sylph-profile", required=True, type=Path)
    parser.add_argument("--vpd-lookup", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    arguments = parser.parse_args()

    sylph_profile = pl.read_csv(arguments.sylph_profile, separator="\t", infer_schema=False, quote_char=None)
    vpd_lookup = pl.read_csv(arguments.vpd_lookup, separator="\t", infer_schema=False, quote_char=None)
    profile_genome_files = set(sylph_profile["Genome_file"].str.split("/").list.last())
    missing_genome_files = sorted(profile_genome_files - set(vpd_lookup["genome_file"]))
    if missing_genome_files:
        sys.stderr.write(f"sylph genome files not in the VPD lookup: {', '.join(missing_genome_files)}\n")
        raise SystemExit(1)

    calls = build_genome_kmer_based_calls(arguments.run_id, arguments.sample_id, sylph_profile, vpd_lookup)
    arguments.out_dir.mkdir(parents=True, exist_ok=True)
    calls.write_csv(arguments.out_dir / f"{arguments.run_id}_{arguments.sample_id}.genome_kmer_based_calls.tsv",
                    separator="\t", null_value="")


if __name__ == "__main__":
    main()
