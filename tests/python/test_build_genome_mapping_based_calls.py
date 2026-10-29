"""Seam A: bin/build_genome_mapping_based_calls.py turns EsViritu's tables into the LabKey call table.

Fixtures are real EsViritu 1.3.5 tables from the ESVIRITU nf-test; expected values are read off
them by hand (e.g. PMMoV 800 reads + ToBRFV 400 reads = the Tobamovirus total of 1,200).
"""

import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

BUILD_SCRIPT = Path(__file__).parents[2] / "bin" / "build_genome_mapping_based_calls.py"
FIXTURES = Path(__file__).parent / "fixtures" / "esviritu"
VPD_LOOKUP = Path(__file__).parents[1] / "data" / "generated" / "references" / "vpd_lookup.tsv"
COLUMNS = [
    "run_id", "sample_id", "vpd_assembly", "call_level", "vpd_species", "vpd_subspecies", "vpd_genus",
    "vpd_family", "vpd_taxids", "subspecies_label", "read_count", "covered_bases", "assembly_length",
    "breadth_fraction", "mean_coverage", "avg_read_identity", "consensus_ref_identity", "rpkmf",
    "normalized_read_ratio",
]


def build(tmp_path: Path, sample_id: str, status: str = "ok", min_read_count: int = 1) -> pl.DataFrame:
    table_arguments = [] if status != "ok" else [
        "--assembly-summary", str(FIXTURES / f"{sample_id}.detected_virus.assembly_summary.tsv"),
        "--info", str(FIXTURES / f"{sample_id}.detected_virus.info.tsv"),
    ]
    completed = subprocess.run(
        [sys.executable, str(BUILD_SCRIPT), "--run-id", "rww_20261029_test", "--sample-id", sample_id,
         "--mapping-status", status, *table_arguments, "--vpd-lookup", str(VPD_LOOKUP),
         "--min-read-count", str(min_read_count), "--species-threshold", "0.9",
         "--subspecies-threshold", "0.95", "--normalization-genus", "Tobamovirus", "--out-dir", str(tmp_path)],
        capture_output=True, text=True,
    )
    assert completed.returncode == 0, completed.stderr
    table = pl.read_csv(tmp_path / f"rww_20261029_test_{sample_id}.genome_mapping_based_calls.tsv",
                        separator="\t", infer_schema=False)
    assert table.columns == COLUMNS
    return table


def row_for(table: pl.DataFrame, assembly: str) -> dict:
    return table.filter(pl.col("vpd_assembly") == assembly).row(0, named=True)


def test_species_level_call_carries_vpd_taxonomy_and_coverage(tmp_path):
    pmmov = row_for(build(tmp_path, "both_callers"), "GCF_000859645.1")

    assert {key: pmmov[key] for key in COLUMNS[:12]} == {
        "run_id": "rww_20261029_test", "sample_id": "both_callers", "vpd_assembly": "GCF_000859645.1",
        "call_level": "species", "vpd_species": "Tobamovirus capsici", "vpd_subspecies": "Pepper mild mottle virus",
        "vpd_genus": "Tobamovirus", "vpd_family": "Virgaviridae", "vpd_taxids": "12239",
        "subspecies_label": "Pepper mild mottle virus", "read_count": "800", "covered_bases": "6351",
    }
    assert pmmov["assembly_length"] == "6357"
    assert float(pmmov["breadth_fraction"]) == pytest.approx(6351 / 6357)
    assert float(pmmov["mean_coverage"]) == pytest.approx(18.876828692779615)
    assert float(pmmov["rpkmf"]) == pytest.approx(33117.243320665984)


def test_normalized_read_ratio_divides_by_the_genus_total(tmp_path):
    table = build(tmp_path, "both_callers")

    assert float(row_for(table, "GCF_000859645.1")["normalized_read_ratio"]) == pytest.approx(800 / 1200)
    assert float(row_for(table, "GCF_001461485.1")["normalized_read_ratio"]) == pytest.approx(400 / 1200)
    assert float(row_for(table, "acc:ON897592.1")["normalized_read_ratio"]) == pytest.approx(956 / 1200)


def test_segmented_assembly_is_one_row_with_assembly_wide_breadth(tmp_path):
    rotavirus = row_for(build(tmp_path, "both_callers"), "set:43d19931109e46219fb3876e13775cee063e1ed0")

    assert float(rotavirus["breadth_fraction"]) == pytest.approx(9997 / 10370)
    # Length-weighted across 11 segments, so it sits between the lowest and highest segment means.
    assert 6.09 < float(rotavirus["mean_coverage"]) < 9.91


def test_below_the_species_threshold_is_a_genus_level_call_with_nothing_below_genus(tmp_path):
    table = build(tmp_path, "mapping_only")

    norovirus = row_for(table, "GCF_001595675.1")
    assert norovirus["call_level"] == "genus"
    assert norovirus["vpd_genus"] == "Norovirus"
    assert norovirus["vpd_species"] is None
    assert norovirus["vpd_subspecies"] is None
    assert norovirus["vpd_taxids"] is None
    assert norovirus["subspecies_label"] is None
    # No Tobamovirus call in this sample, so there is no denominator.
    assert norovirus["normalized_read_ratio"] is None


def test_calls_under_the_minimum_read_count_are_dropped_before_normalising(tmp_path):
    table = build(tmp_path, "both_callers", min_read_count=700)

    assert sorted(table["vpd_assembly"]) == ["GCF_000859645.1", "acc:ON897592.1"]
    assert float(row_for(table, "GCF_000859645.1")["normalized_read_ratio"]) == pytest.approx(1.0)


@pytest.mark.parametrize("status", ["no_detections", "failed"])
def test_no_tables_from_esviritu_gives_a_header_only_table(tmp_path, status):
    assert build(tmp_path, "kmer_only", status=status).height == 0
