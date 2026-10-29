#!/usr/bin/env python3
"""Build one sample's genome_mapping_based_calls table from EsViritu's output tables.

A genome mapping-based call is one EsViritu Assembly with read_count >= --min-read-count.

call_level: EsViritu relabels low-identity calls to "s__unclassified_<genus>" only in its
tax_profile.tsv (esv_funcs.tax_profile); detected_virus.assembly_summary.tsv keeps the species.
This script applies the same rule to the assembly table: effective identity =
consensus_ref_identity, or avg_read_identity when that is missing; below --species-threshold
the call is genus level and every column below genus (vpd_species, vpd_subspecies, vpd_taxids,
subspecies_label) is null. subspecies_label is EsViritu's subspecies label, null below
--subspecies-threshold (EsViritu would label it unclassified).

Taxonomy comes from the VPD lookup (shared with genome_kmer_based_calls), joined on Assembly.
mean_coverage is length-weighted over the Assembly's records from detected_virus.info.tsv, and
breadth_fraction is covered_bases / Asm_length.

normalized_read_ratio = read_count / the sample's total read_count over calls whose vpd_genus is
--normalization-genus (any call_level); null when the sample has no call in that genus.

With --mapping-status no_detections or failed, EsViritu wrote no tables and the output is a
header-only table; the status itself goes to sample_status.

Output: <out-dir>/<run_id>_<sample_id>.genome_mapping_based_calls.tsv
"""

import argparse
from pathlib import Path

import polars as pl

COLUMN_TYPES = {
    "run_id": pl.String, "sample_id": pl.String, "vpd_assembly": pl.String, "call_level": pl.String,
    "vpd_species": pl.String, "vpd_subspecies": pl.String, "vpd_genus": pl.String, "vpd_family": pl.String,
    "vpd_taxids": pl.String, "subspecies_label": pl.String, "read_count": pl.Int64, "covered_bases": pl.Int64,
    "assembly_length": pl.Int64, "breadth_fraction": pl.Float64, "mean_coverage": pl.Float64,
    "avg_read_identity": pl.Float64, "consensus_ref_identity": pl.Float64, "rpkmf": pl.Float64,
    "normalized_read_ratio": pl.Float64,
}


def read_tsv(path: Path) -> pl.DataFrame:
    return pl.read_csv(path, separator="\t", infer_schema=False, quote_char=None)


def assembly_mean_coverage(info_table: pl.DataFrame) -> pl.DataFrame:
    """Length-weighted mean coverage per Assembly over its records (segments)."""
    return info_table.group_by("Assembly").agg(
        ((pl.col("mean_coverage").cast(pl.Float64) * pl.col("Length").cast(pl.Float64)).sum()
         / pl.col("Asm_length").cast(pl.Float64).first()).alias("mean_coverage"))


def build_genome_mapping_based_calls(run_id: str, sample_id: str, assembly_summary: pl.DataFrame,
                                     info_table: pl.DataFrame, vpd_lookup: pl.DataFrame, min_read_count: int,
                                     species_threshold: float, subspecies_threshold: float,
                                     normalization_genus: str) -> pl.DataFrame:
    effective_identity = pl.coalesce(pl.col("consensus_ref_identity"), pl.col("avg_read_identity"))
    is_genus_level = effective_identity < species_threshold
    calls = (
        assembly_summary
        .with_columns(
            pl.col("read_count").cast(pl.Int64),
            pl.col("covered_bases").cast(pl.Int64),
            pl.col("Asm_length").cast(pl.Int64).alias("assembly_length"),
            pl.col("avg_read_identity").cast(pl.Float64),
            pl.col("consensus_ref_identity").cast(pl.Float64, strict=False),
            pl.col("RPKMF").cast(pl.Float64).alias("rpkmf"),
        )
        .filter(pl.col("read_count") >= min_read_count)
        .join(assembly_mean_coverage(info_table), on="Assembly", how="left")
        .join(vpd_lookup, left_on="Assembly", right_on="vpd_assembly", how="left")
        .with_columns(
            pl.col("Assembly").alias("vpd_assembly"),
            pl.when(is_genus_level).then(pl.lit("genus")).otherwise(pl.lit("species")).alias("call_level"),
            pl.when(is_genus_level).then(None).otherwise(pl.col("vpd_species")).alias("vpd_species"),
            pl.when(is_genus_level).then(None).otherwise(pl.col("vpd_subspecies")).alias("vpd_subspecies"),
            pl.when(is_genus_level).then(None).otherwise(pl.col("vpd_taxids")).alias("vpd_taxids"),
            pl.when(is_genus_level | (effective_identity < subspecies_threshold)).then(None)
            .otherwise(pl.col("subspecies").str.replace(r"^t__", "")).alias("subspecies_label"),
            (pl.col("covered_bases") / pl.col("assembly_length")).alias("breadth_fraction"),
        )
    )
    genus_read_total = calls.filter(pl.col("vpd_genus") == normalization_genus)["read_count"].sum()
    return calls.with_columns(
        pl.lit(run_id).alias("run_id"),
        pl.lit(sample_id).alias("sample_id"),
        (pl.col("read_count") / genus_read_total if genus_read_total else pl.lit(None, pl.Float64))
        .alias("normalized_read_ratio"),
    ).select([pl.col(name).cast(dtype) for name, dtype in COLUMN_TYPES.items()])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--mapping-status", required=True, choices=["ok", "no_detections", "failed"])
    parser.add_argument("--assembly-summary", type=Path, help="<sample>.detected_virus.assembly_summary.tsv")
    parser.add_argument("--info", type=Path, help="<sample>.detected_virus.info.tsv")
    parser.add_argument("--vpd-lookup", required=True, type=Path)
    parser.add_argument("--min-read-count", required=True, type=int)
    parser.add_argument("--species-threshold", required=True, type=float)
    parser.add_argument("--subspecies-threshold", required=True, type=float)
    parser.add_argument("--normalization-genus", required=True)
    parser.add_argument("--out-dir", required=True, type=Path)
    arguments = parser.parse_args()

    if arguments.mapping_status == "ok":
        calls = build_genome_mapping_based_calls(
            arguments.run_id, arguments.sample_id, read_tsv(arguments.assembly_summary), read_tsv(arguments.info),
            read_tsv(arguments.vpd_lookup), arguments.min_read_count, arguments.species_threshold,
            arguments.subspecies_threshold, arguments.normalization_genus)
    else:
        calls = pl.DataFrame(schema=COLUMN_TYPES)
    arguments.out_dir.mkdir(parents=True, exist_ok=True)
    calls.write_csv(arguments.out_dir / f"{arguments.run_id}_{arguments.sample_id}.genome_mapping_based_calls.tsv",
                    separator="\t", null_value="")


if __name__ == "__main__":
    main()
