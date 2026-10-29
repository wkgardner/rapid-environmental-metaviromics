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
