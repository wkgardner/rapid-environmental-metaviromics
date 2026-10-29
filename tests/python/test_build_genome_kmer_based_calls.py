"""Seam A: bin/build_genome_kmer_based_calls.py turns a sylph profile into the LabKey call table.

Fixtures are real sylph 1.0.0 profiles from the SYLPH_PROFILE nf-test; the lookup is the
committed test reference built from the VPD subset.
"""

import subprocess
import sys
from pathlib import Path

import polars as pl

BUILD_SCRIPT = Path(__file__).parents[2] / "bin" / "build_genome_kmer_based_calls.py"
FIXTURES = Path(__file__).parent / "fixtures" / "sylph"
VPD_LOOKUP = Path(__file__).parents[1] / "data" / "generated" / "references" / "vpd_lookup.tsv"
COLUMNS = [
    "run_id", "sample_id", "vpd_assembly", "vpd_species", "vpd_subspecies", "vpd_genus", "vpd_family",
    "vpd_taxids", "taxonomic_abundance", "sequence_abundance", "adjusted_ani", "naive_ani",
    "effective_coverage", "containment_matching_kmers", "containment_total_kmers", "kmers_reassigned",
]


def build(tmp_path: Path, sample_id: str, profile: Path, lookup: Path = VPD_LOOKUP) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(BUILD_SCRIPT), "--run-id", "rww_20261029_test", "--sample-id", sample_id,
         "--sylph-profile", str(profile), "--vpd-lookup", str(lookup), "--out-dir", str(tmp_path)],
        capture_output=True, text=True,
    )


def read_table(tmp_path: Path, sample_id: str) -> pl.DataFrame:
    return pl.read_csv(tmp_path / f"rww_20261029_test_{sample_id}.genome_kmer_based_calls.tsv",
                       separator="\t", infer_schema=False)


def test_one_row_per_profile_row_joined_to_vpd_by_genome_file(tmp_path):
    completed = build(tmp_path, "both_callers", FIXTURES / "both_callers.sylph_profile.tsv")
    assert completed.returncode == 0, completed.stderr

    table = read_table(tmp_path, "both_callers")
    assert table.columns == COLUMNS
    assert sorted(table["vpd_assembly"]) == sorted([
        "GCF_000859645.1", "GCF_001461485.1", "acc:ON897592.1",
        "set:43d19931109e46219fb3876e13775cee063e1ed0", "PSVA/S2-2__Sadwavirus_alphananas",
    ])
    pmmov = table.filter(pl.col("vpd_assembly") == "GCF_000859645.1").row(0, named=True)
    assert pmmov == {
        "run_id": "rww_20261029_test", "sample_id": "both_callers", "vpd_assembly": "GCF_000859645.1",
        "vpd_species": "Tobamovirus capsici", "vpd_subspecies": "Pepper mild mottle virus",
        "vpd_genus": "Tobamovirus", "vpd_family": "Virgaviridae", "vpd_taxids": "12239",
        "taxonomic_abundance": "35.4178", "sequence_abundance": "22.7888", "adjusted_ani": "100.0",
        "naive_ani": "100.0", "effective_coverage": "15.0",
        "containment_matching_kmers": "81", "containment_total_kmers": "81", "kmers_reassigned": "0",
    }


def test_a_profile_with_no_rows_gives_a_header_only_table(tmp_path):
    completed = build(tmp_path, "mapping_only", FIXTURES / "mapping_only.sylph_profile.tsv")
    assert completed.returncode == 0, completed.stderr

    table = read_table(tmp_path, "mapping_only")
    assert table.columns == COLUMNS
    assert table.height == 0


def test_a_genome_file_missing_from_the_lookup_stops_the_build(tmp_path):
    partial_lookup = tmp_path / "partial_lookup.tsv"
    pl.read_csv(VPD_LOOKUP, separator="\t", infer_schema=False).filter(
        pl.col("vpd_assembly") != "GCF_000859645.1").write_csv(partial_lookup, separator="\t")

    completed = build(tmp_path, "both_callers", FIXTURES / "both_callers.sylph_profile.tsv", partial_lookup)
    assert completed.returncode == 1
    assert "GCF_000859645.1.fna" in completed.stderr
