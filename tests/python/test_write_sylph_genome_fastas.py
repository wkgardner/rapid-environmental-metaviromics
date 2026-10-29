"""Seam A: bin/write_sylph_genome_fastas.py writes one N-free FASTA per VPD Assembly for sylph."""

import subprocess
import sys
from pathlib import Path

from Bio import SeqIO

WRITE_SCRIPT = Path(__file__).parents[2] / "bin" / "write_sylph_genome_fastas.py"
LOOKUP_SCRIPT = Path(__file__).parents[2] / "bin" / "build_vpd_lookup.py"
VPD_SUBSET = Path(__file__).parents[1] / "data" / "vpd_subset"


def write_genomes(tmp_path: Path) -> Path:
    lookup_path = tmp_path / "vpd_lookup.tsv"
    subprocess.run([sys.executable, str(LOOKUP_SCRIPT), "--vpd-metadata",
                    str(VPD_SUBSET / "virus_pathogen_database.all_metadata.tsv"), "--out", str(lookup_path)], check=True)
    completed = subprocess.run(
        [sys.executable, str(WRITE_SCRIPT),
         "--vpd-fasta", str(VPD_SUBSET / "virus_pathogen_database.fna"),
         "--vpd-metadata", str(VPD_SUBSET / "virus_pathogen_database.all_metadata.tsv"),
         "--vpd-lookup", str(lookup_path), "--min-fragment-length", "31",
         "--out-dir", str(tmp_path / "genomes"), "--genome-list", str(tmp_path / "genome_list.txt")],
        capture_output=True, text=True,
    )
    assert completed.returncode == 0, completed.stderr
    return tmp_path / "genomes"


def test_one_fasta_per_assembly_listed_in_the_genome_list(tmp_path):
    genomes = write_genomes(tmp_path)

    listed = (tmp_path / "genome_list.txt").read_text().split()
    assert len(listed) == 7
    assert sorted(Path(path).name for path in listed) == sorted(path.name for path in genomes.iterdir())
    assert (genomes / "PSVA_S2-2__Sadwavirus_alphananas.fna").exists()


def test_segments_stay_together_in_their_assembly_file(tmp_path):
    genomes = write_genomes(tmp_path)

    rotavirus = list(SeqIO.parse(genomes / "set_43d19931109e46219fb3876e13775cee063e1ed0.fna", "fasta"))
    assert len(rotavirus) == 11


def test_records_are_cut_at_n_runs_and_short_fragments_dropped(tmp_path):
    genomes = write_genomes(tmp_path)

    # MF678311.1 is 7,183 bp with one N at position 4,747.
    enterovirus = list(SeqIO.parse(genomes / "set_MF678311.fna", "fasta"))
    assert [(record.id, len(record.seq)) for record in enterovirus] == [("MF678311.1__f0", 4747), ("MF678311.1__f1", 2435)]
    for path in genomes.iterdir():
        for record in SeqIO.parse(path, "fasta"):
            assert "N" not in str(record.seq).upper()
            assert len(record.seq) >= 31
