#!/usr/bin/env python3
"""Write a seeded synthetic dataset for the pipeline's tests from a YAML spec.

The spec names samples and, for each, the read sources mixed into it. A source is
a VPD Assembly from the committed subset, or the synthetic "host", "background" or
"low_complexity" (dinucleotide repeat) genome. Outputs, under --out-dir:

    samplesheet.csv                      NVD-format samplesheet, one row per sample
    reads/<sample_id>_R1/_R2.fastq.gz    paired reads (gzip written with mtime 0)
    truth/<sample_id>.reads.tsv          true origin of every pair
    references/host.fna                  synthetic host genome, for a depletion index

The same spec and seed always give byte-identical files.
"""

import argparse
import gzip
import random
import re
from pathlib import Path

import polars as pl
import yaml
from Bio import SeqIO

NUCLEOTIDES = "ACGT"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
DEFAULT_READ_LENGTH = 150
DEFAULT_MEAN_FRAGMENT_LENGTH = 300
FRAGMENT_LENGTH_SPREAD = 50
HIGH_QUALITY_CHARACTER = "I"  # Phred 40
LOW_QUALITY_CHARACTER = "#"  # Phred 2, so a quality filter removes the read
# Illumina TruSeq adapters: what a read runs into when its insert is shorter than the read.
TRUSEQ_ADAPTER_READ1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
TRUSEQ_ADAPTER_READ2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
SHORTEST_ADAPTER_INSERT = 30
NVD_SAMPLESHEET_COLUMNS = ["sample_id", "srr", "platform", "fastq1", "fastq2", "fastq1_glob", "fastq2_glob"]


def sanitise_label(name: str) -> str:
    """Make a genome name safe for read names, e.g. 'set:MF678311' -> 'set_MF678311'."""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def random_genome(random_generator: random.Random, length: int) -> str:
    return "".join(random_generator.choice(NUCLEOTIDES) for _ in range(length))


def load_vpd_assembly_records(vpd_fasta: Path, vpd_metadata: Path) -> dict[str, list[tuple[str, str]]]:
    """Map each VPD Assembly to its (accession, sequence) records, e.g. all 11 rotavirus segments."""
    sequence_by_accession = {record.id: str(record.seq).upper() for record in SeqIO.parse(vpd_fasta, "fasta")}
    metadata = pl.read_csv(vpd_metadata, separator="\t", infer_schema=False, quote_char=None)
    records_by_assembly: dict[str, list[tuple[str, str]]] = {}
    for row in metadata.iter_rows(named=True):
        records_by_assembly.setdefault(row["Assembly"], []).append((row["Accession"], sequence_by_accession[row["Accession"]]))
    return records_by_assembly


def read_through_adapter(insert: str, adapter: str, read_length: int) -> str:
    """A read longer than its insert reads the insert, then the adapter, then a poly-A tail."""
    return (insert + adapter + "A" * read_length)[:read_length]


def low_complexity_genome(random_generator: random.Random, length: int) -> str:
    """A dinucleotide repeat such as ACACAC..., which an entropy filter should remove."""
    first, second = random_generator.sample(NUCLEOTIDES, 2)
    return ((first + second) * length)[:length]


def mutate(sequence: str, divergence: float, random_generator: random.Random) -> str:
    """Substitute each base with probability `divergence` (always to a different base)."""
    if divergence == 0:
        return sequence
    bases = list(sequence)
    for position, base in enumerate(bases):
        if random_generator.random() < divergence:
            bases[position] = random_generator.choice([other for other in NUCLEOTIDES if other != base])
    return "".join(bases)


def simulate_source_pairs(source: dict, records: list[tuple[str, str]], read_name_prefix: str,
                          random_generator: random.Random) -> list[dict]:
    """Simulate one source's read pairs: unique pairs first, then exact copies for duplicates."""
    read_length = source.get("read_length", DEFAULT_READ_LENGTH)
    quality_character = LOW_QUALITY_CHARACTER if source.get("low_quality", False) else HIGH_QUALITY_CHARACTER
    total_pairs = source["read_pairs"]
    duplicate_pairs = round(total_pairs * source.get("duplicate_fraction", 0.0))
    record_lengths = [len(sequence) for _, sequence in records]

    pairs = []
    adapter_pairs = round((total_pairs - duplicate_pairs) * source.get("adapter_fraction", 0.0))
    # A diverged source is one strain: mutate each record once, then sample reads from it exactly,
    # so every read carries the same mutations (independent per-read errors would let k-mer tools
    # recover the reference k-mers at high coverage).
    divergence = source.get("divergence", 0.0)
    records = [(accession, mutate(sequence, divergence, random_generator)) for accession, sequence in records]
    for pair_index in range(total_pairs - duplicate_pairs):
        accession, sequence = random_generator.choices(records, weights=record_lengths)[0]
        reads_into_adapter = pair_index < adapter_pairs
        if reads_into_adapter:
            fragment_length = random_generator.randrange(SHORTEST_ADAPTER_INSERT, read_length - SHORTEST_ADAPTER_INSERT)
        else:
            fragment_length = min(len(sequence), max(read_length, int(random_generator.gauss(
                DEFAULT_MEAN_FRAGMENT_LENGTH, FRAGMENT_LENGTH_SPREAD))))
        fragment_start = random_generator.randrange(0, len(sequence) - fragment_length + 1)
        fragment = sequence[fragment_start: fragment_start + fragment_length]
        if reads_into_adapter:
            read1 = read_through_adapter(fragment, TRUSEQ_ADAPTER_READ1, read_length)
            read2 = read_through_adapter(reverse_complement(fragment), TRUSEQ_ADAPTER_READ2, read_length)
        else:
            read1 = fragment[:read_length]
            read2 = reverse_complement(fragment)[:read_length]
        pairs.append({
            "source_record": accession,
            "read1_start": fragment_start,
            "fragment_length": fragment_length,
            "read1": read1,
            "read2": read2,
            "quality": quality_character * read_length,
            "is_duplicate": False,
        })
    unique_pairs = list(pairs)
    for _ in range(duplicate_pairs):
        pairs.append(dict(random_generator.choice(unique_pairs), is_duplicate=True))

    for index, pair in enumerate(pairs):
        pair["read_name"] = f"{read_name_prefix}_{index}"
    return pairs


def write_gzip_fastq(path: Path, records: list[tuple[str, str, str]]) -> None:
    """Write FASTQ through gzip with mtime 0 so reruns are byte-identical."""
    text = "".join(f"@{name}\n{sequence}\n+\n{quality}\n" for name, sequence, quality in records)
    with open(path, "wb") as raw_handle, gzip.GzipFile(fileobj=raw_handle, mode="wb", mtime=0) as handle:
        handle.write(text.encode())


def make_test_data(spec_path: Path, out_dir: Path) -> None:
    spec = yaml.safe_load(Path(spec_path).read_text())
    random_generator = random.Random(spec["seed"])
    out_dir = Path(out_dir).resolve()
    for subdirectory in ("reads", "truth", "references"):
        (out_dir / subdirectory).mkdir(parents=True, exist_ok=True)

    genomes = load_vpd_assembly_records(Path(spec["vpd_fasta"]), Path(spec["vpd_metadata"]))
    genomes["host"] = [("synthetic_host", random_genome(random_generator, spec["host_genome_length"]))]
    genomes["low_complexity"] = [("synthetic_low_complexity", low_complexity_genome(random_generator, 5_000))]
    genomes["background"] = [("synthetic_background", random_genome(random_generator, spec["background_genome_length"]))]
    (out_dir / "references" / "host.fna").write_text(f">synthetic_host\n{genomes['host'][0][1]}\n")

    samplesheet_rows = []
    for sample in spec["samples"]:
        sample_id = sample["sample_id"]
        pairs = []
        for source in sample["sources"]:
            prefix = f"{sample_id}_{sanitise_label(source['genome'])}"
            pairs.extend(simulate_source_pairs(source, genomes[source["genome"]], prefix, random_generator))

        read1_path = out_dir / "reads" / f"{sample_id}_R1.fastq.gz"
        read2_path = out_dir / "reads" / f"{sample_id}_R2.fastq.gz"
        write_gzip_fastq(read1_path, [(f"{pair['read_name']}/1", pair["read1"], pair["quality"]) for pair in pairs])
        write_gzip_fastq(read2_path, [(f"{pair['read_name']}/2", pair["read2"], pair["quality"]) for pair in pairs])
        pl.DataFrame(
            [{key: pair[key] for key in ("read_name", "source_record", "read1_start", "fragment_length", "is_duplicate")}
             for pair in pairs],
            schema={"read_name": pl.String, "source_record": pl.String, "read1_start": pl.Int64,
                    "fragment_length": pl.Int64, "is_duplicate": pl.Boolean},
        ).write_csv(out_dir / "truth" / f"{sample_id}.reads.tsv", separator="\t")
        samplesheet_rows.append({"sample_id": sample_id, "srr": "", "platform": "illumina",
                                 "fastq1": str(read1_path), "fastq2": str(read2_path),
                                 "fastq1_glob": "", "fastq2_glob": ""})

    pl.DataFrame(samplesheet_rows, schema={column: pl.String for column in NVD_SAMPLESHEET_COLUMNS}).write_csv(
        out_dir / "samplesheet.csv")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--spec", required=True, type=Path, help="YAML spec of samples and read sources")
    parser.add_argument("--out-dir", required=True, type=Path, help="directory to write the dataset into")
    arguments = parser.parse_args()
    make_test_data(arguments.spec, arguments.out_dir)


if __name__ == "__main__":
    main()
