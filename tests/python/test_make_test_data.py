"""Seam A: bin/make_test_data.py writes a seeded synthetic dataset from a YAML spec."""

import gzip
import hashlib
from pathlib import Path

import polars as pl
import pytest
import yaml
from Bio import SeqIO

from make_test_data import make_test_data

VPD_SUBSET_DIRECTORY = Path(__file__).parents[1] / "data" / "vpd_subset"


def write_spec(tmp_path: Path, samples: list[dict], seed: int = 7) -> Path:
    """Write a small generator spec that uses the committed VPD subset."""
    spec = {
        "seed": seed,
        "vpd_fasta": str(VPD_SUBSET_DIRECTORY / "virus_pathogen_database.fna"),
        "vpd_metadata": str(VPD_SUBSET_DIRECTORY / "virus_pathogen_database.all_metadata.tsv"),
        "host_genome_length": 20_000,
        "background_genome_length": 10_000,
        "samples": samples,
    }
    spec_path = tmp_path / "spec.yaml"
    spec_path.write_text(yaml.safe_dump(spec))
    return spec_path


def read_fastq_records(fastq_path: Path) -> list[tuple[str, str, str]]:
    """Return (name, sequence, quality) for every record in a gzipped FASTQ."""
    with gzip.open(fastq_path, "rt") as handle:
        lines = handle.read().splitlines()
    return [(lines[i][1:], lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4)]


def file_digest(path: Path) -> str:
    return hashlib.sha256(gzip.open(path).read()).hexdigest()


TWO_SAMPLES = [
    {"sample_id": "both_callers", "sources": [
        {"genome": "GCF_000859645.1", "read_pairs": 50},
        {"genome": "host", "read_pairs": 30},
    ]},
    {"sample_id": "empty_sample", "sources": []},
]


def test_samplesheet_has_one_nvd_row_per_sample(tmp_path):
    out_dir = tmp_path / "out"
    make_test_data(write_spec(tmp_path, TWO_SAMPLES), out_dir)

    samplesheet = pl.read_csv(out_dir / "samplesheet.csv", infer_schema=False)
    assert samplesheet.columns == ["sample_id", "srr", "platform", "fastq1", "fastq2", "fastq1_glob", "fastq2_glob"]
    assert samplesheet["sample_id"].to_list() == ["both_callers", "empty_sample"]
    assert samplesheet["platform"].to_list() == ["illumina", "illumina"]
    assert Path(samplesheet["fastq1"][0]).name == "both_callers_R1.fastq.gz"


def test_each_sample_gets_the_requested_number_of_pairs(tmp_path):
    out_dir = tmp_path / "out"
    make_test_data(write_spec(tmp_path, TWO_SAMPLES), out_dir)

    read1 = read_fastq_records(out_dir / "reads" / "both_callers_R1.fastq.gz")
    read2 = read_fastq_records(out_dir / "reads" / "both_callers_R2.fastq.gz")
    assert len(read1) == len(read2) == 80
    assert [name.split("/")[0] for name, _, _ in read1] == [name.split("/")[0] for name, _, _ in read2]
    assert read_fastq_records(out_dir / "reads" / "empty_sample_R1.fastq.gz") == []


def test_same_seed_gives_identical_reads_and_a_new_seed_does_not(tmp_path):
    make_test_data(write_spec(tmp_path, TWO_SAMPLES, seed=7), tmp_path / "first")
    make_test_data(write_spec(tmp_path, TWO_SAMPLES, seed=7), tmp_path / "second")
    make_test_data(write_spec(tmp_path, TWO_SAMPLES, seed=8), tmp_path / "third")

    read1_name = Path("reads") / "both_callers_R1.fastq.gz"
    assert file_digest(tmp_path / "first" / read1_name) == file_digest(tmp_path / "second" / read1_name)
    assert file_digest(tmp_path / "first" / read1_name) != file_digest(tmp_path / "third" / read1_name)


def test_read_length_and_quality_follow_the_source(tmp_path):
    samples = [{"sample_id": "short_reads", "sources": [
        {"genome": "set:MF678311", "read_pairs": 20, "read_length": 80},
        {"genome": "background", "read_pairs": 10, "low_quality": True},
    ]}]
    out_dir = tmp_path / "out"
    make_test_data(write_spec(tmp_path, samples), out_dir)

    records = read_fastq_records(out_dir / "reads" / "short_reads_R1.fastq.gz")
    enterovirus_reads = [record for record in records if "set_MF678311" in record[0]]
    background_reads = [record for record in records if "background" in record[0]]
    assert {len(sequence) for _, sequence, _ in enterovirus_reads} == {80}
    assert {quality for _, _, quality in enterovirus_reads} == {"I" * 80}
    assert {character for _, _, quality in background_reads for character in quality} == {"#"}


@pytest.mark.parametrize("divergence", [0.0, 0.12])
def test_divergence_sets_the_mismatch_rate_against_the_reference(tmp_path, divergence):
    samples = [{"sample_id": "diverged", "sources": [
        {"genome": "GCF_001595675.1", "read_pairs": 200, "divergence": divergence},
    ]}]
    out_dir = tmp_path / "out"
    make_test_data(write_spec(tmp_path, samples), out_dir)

    # The generator writes each read's true origin, so the reference base under every read is known.
    truth = pl.read_csv(out_dir / "truth" / "diverged.reads.tsv", separator="\t")
    reference_sequences = {
        record.id: str(record.seq)
        for record in SeqIO.parse(VPD_SUBSET_DIRECTORY / "virus_pathogen_database.fna", "fasta")
    }
    records = {name.split("/")[0]: sequence for name, sequence, _ in read_fastq_records(out_dir / "reads" / "diverged_R1.fastq.gz")}

    mismatches = compared = 0
    for row in truth.iter_rows(named=True):
        read_sequence = records[row["read_name"]]
        source_sequence = reference_sequences[row["source_record"]]
        reference_window = source_sequence[row["read1_start"]: row["read1_start"] + len(read_sequence)]
        compared += len(read_sequence)
        mismatches += sum(a != b for a, b in zip(read_sequence, reference_window))
    assert mismatches / compared == pytest.approx(divergence, abs=0.02)


def test_duplicate_fraction_writes_exact_copies_of_pairs(tmp_path):
    samples = [{"sample_id": "duplicated", "sources": [
        {"genome": "acc:ON897592.1", "read_pairs": 100, "duplicate_fraction": 0.3},
    ]}]
    out_dir = tmp_path / "out"
    make_test_data(write_spec(tmp_path, samples), out_dir)

    read1 = read_fastq_records(out_dir / "reads" / "duplicated_R1.fastq.gz")
    read2 = read_fastq_records(out_dir / "reads" / "duplicated_R2.fastq.gz")
    pair_sequences = [(first[1], second[1]) for first, second in zip(read1, read2)]
    assert len(pair_sequences) == 100
    assert len(pair_sequences) - len(set(pair_sequences)) == 30


def test_host_genome_is_written_for_building_a_depletion_index(tmp_path):
    out_dir = tmp_path / "out"
    make_test_data(write_spec(tmp_path, TWO_SAMPLES), out_dir)

    host_fasta = (out_dir / "references" / "host.fna").read_text().splitlines()
    assert host_fasta[0] == ">synthetic_host"
    assert len("".join(host_fasta[1:])) == 20_000


TRUSEQ_ADAPTER_READ1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"


def test_adapter_fraction_gives_short_inserts_that_read_into_the_adapter(tmp_path):
    samples = [{"sample_id": "adapters", "sources": [
        {"genome": "acc:ON897592.1", "read_pairs": 40, "adapter_fraction": 0.5},
    ]}]
    out_dir = tmp_path / "out"
    make_test_data(write_spec(tmp_path, samples), out_dir)

    truth = pl.read_csv(out_dir / "truth" / "adapters.reads.tsv", separator="\t")
    read1 = {name.split("/")[0]: sequence for name, sequence, _ in read_fastq_records(out_dir / "reads" / "adapters_R1.fastq.gz")}
    short_inserts = truth.filter(pl.col("fragment_length") < 150)
    assert short_inserts.height == 20
    for row in short_inserts.iter_rows(named=True):
        sequence = read1[row["read_name"]]
        assert len(sequence) == 150
        assert sequence[row["fragment_length"]:].startswith(TRUSEQ_ADAPTER_READ1[:20])


def test_low_complexity_source_writes_dinucleotide_repeats(tmp_path):
    samples = [{"sample_id": "repeats", "sources": [{"genome": "low_complexity", "read_pairs": 10}]}]
    out_dir = tmp_path / "out"
    make_test_data(write_spec(tmp_path, samples), out_dir)

    for _, sequence, _ in read_fastq_records(out_dir / "reads" / "repeats_R1.fastq.gz"):
        assert len(set(sequence)) <= 2


def test_a_diverged_source_is_one_strain_so_overlapping_reads_agree(tmp_path):
    samples = [{"sample_id": "strain", "sources": [
        {"genome": "GCF_001595675.1", "read_pairs": 300, "divergence": 0.12},
    ]}]
    out_dir = tmp_path / "out"
    make_test_data(write_spec(tmp_path, samples), out_dir)

    truth = pl.read_csv(out_dir / "truth" / "strain.reads.tsv", separator="\t")
    read1 = {name.split("/")[0]: sequence for name, sequence, _ in read_fastq_records(out_dir / "reads" / "strain_R1.fastq.gz")}
    bases_seen_at_position: dict[int, set[str]] = {}
    for row in truth.iter_rows(named=True):
        for offset, base in enumerate(read1[row["read_name"]]):
            bases_seen_at_position.setdefault(row["read1_start"] + offset, set()).add(base)
    covered_more_than_once = [bases for bases in bases_seen_at_position.values() if bases]
    assert all(len(bases) == 1 for bases in covered_more_than_once)
