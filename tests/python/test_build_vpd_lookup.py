"""Seam A: bin/build_vpd_lookup.py writes one row per VPD Assembly, shared by both callers."""

import subprocess
import sys
from pathlib import Path

import polars as pl

BUILD_SCRIPT = Path(__file__).parents[2] / "bin" / "build_vpd_lookup.py"
VPD_METADATA = Path(__file__).parents[1] / "data" / "vpd_subset" / "virus_pathogen_database.all_metadata.tsv"


def build_lookup(tmp_path: Path, metadata: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(BUILD_SCRIPT), "--vpd-metadata", str(metadata), "--out", str(tmp_path / "vpd_lookup.tsv")],
        capture_output=True, text=True,
    )


def test_one_row_per_assembly_with_rank_prefixes_stripped(tmp_path):
    completed = build_lookup(tmp_path, VPD_METADATA)
    assert completed.returncode == 0, completed.stderr

    lookup = pl.read_csv(tmp_path / "vpd_lookup.tsv", separator="\t", infer_schema=False)
    assert lookup.columns == ["vpd_assembly", "genome_file", "vpd_species", "vpd_subspecies", "vpd_genus", "vpd_family", "vpd_taxids"]
    assert lookup.height == 7
    pmmov = lookup.filter(pl.col("vpd_assembly") == "GCF_000859645.1").row(0, named=True)
    assert pmmov == {
        "vpd_assembly": "GCF_000859645.1", "genome_file": "GCF_000859645.1.fna",
        "vpd_species": "Tobamovirus capsici", "vpd_subspecies": "Pepper mild mottle virus",
        "vpd_genus": "Tobamovirus", "vpd_family": "Virgaviridae", "vpd_taxids": "12239",
    }


def test_genome_files_are_sanitised_and_segment_labels_are_joined(tmp_path):
    build_lookup(tmp_path, VPD_METADATA)
    lookup = pl.read_csv(tmp_path / "vpd_lookup.tsv", separator="\t", infer_schema=False)
    genome_file_by_assembly = dict(zip(lookup["vpd_assembly"], lookup["genome_file"]))

    assert genome_file_by_assembly["PSVA/S2-2__Sadwavirus_alphananas"] == "PSVA_S2-2__Sadwavirus_alphananas.fna"
    assert genome_file_by_assembly["acc:ON897592.1"] == "acc_ON897592.1.fna"
    rotavirus = lookup.filter(pl.col("vpd_genus") == "Rotavirus").row(0, named=True)
    assert rotavirus["vpd_subspecies"] == "A1;C1;E1;G11;H1;I12;M1;N1;P[25];R1;T1"
    assert rotavirus["vpd_taxids"] == "28875"


def test_two_assemblies_that_sanitise_to_the_same_file_name_are_rejected(tmp_path):
    metadata = pl.read_csv(VPD_METADATA, separator="\t", infer_schema=False, quote_char=None)
    clashing = pl.concat([
        metadata.head(1).with_columns(pl.lit("set:AB1").alias("Assembly")),
        metadata.head(1).with_columns(pl.lit("set/AB1").alias("Assembly"), pl.lit("XX000001.1").alias("Accession")),
    ])
    clashing_path = tmp_path / "clashing.tsv"
    clashing.write_csv(clashing_path, separator="\t", quote_style="never")

    completed = build_lookup(tmp_path, clashing_path)
    assert completed.returncode == 1
    assert "set_AB1.fna" in completed.stderr
