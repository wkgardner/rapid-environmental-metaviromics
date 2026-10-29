#!/usr/bin/env bash
# Build the tiny real indexes the nf-tests run against, from tests/data/generated/.
# Called by `pixi run make-test-data` after bin/make_test_data.py.
set -euo pipefail

generated_directory="tests/data/generated"
references_directory="${generated_directory}/references"

# Host depletion index: the synthetic host at the pipeline's k31/w15.
deacon index build -k 31 -w 15 --threads 2 --quiet \
    --output "${references_directory}/host.k31w15.idx" \
    "${references_directory}/host.fna"

# VPD lookup and the genome k-mer-based (sylph) database at the spec's k31 / c50.
# Built from inside references/ so the database records relative genome paths.
vpd_subset_directory="$(pwd)/tests/data/vpd_subset"
(
    cd "${references_directory}"
    rm -rf sylph_genomes
    build_vpd_lookup.py --vpd-metadata "${vpd_subset_directory}/virus_pathogen_database.all_metadata.tsv" \
        --out vpd_lookup.tsv
    write_sylph_genome_fastas.py \
        --vpd-fasta "${vpd_subset_directory}/virus_pathogen_database.fna" \
        --vpd-metadata "${vpd_subset_directory}/virus_pathogen_database.all_metadata.tsv" \
        --vpd-lookup vpd_lookup.tsv --min-fragment-length 31 \
        --out-dir sylph_genomes --genome-list sylph_genome_list.txt
    sylph sketch --gl sylph_genome_list.txt -k 31 -c 50 -t 2 -o vpd
    rm -rf sylph_genomes sylph_genome_list.txt
)
