#!/usr/bin/env bash
# Build a k31/w15 Deacon index of the vector-masked VPD (the same FASTA EsViritu searches), e.g.
# for extracting or screening viral reads. Run from the repo root:
#   pixi run bash scripts/build_vpd_deacon_index.sh <vpd_fasta> <output_index>
set -euo pipefail

vpd_fasta="${1:?usage: build_vpd_deacon_index.sh <vpd_fasta> <output_index>}"
output_index="${2:?usage: build_vpd_deacon_index.sh <vpd_fasta> <output_index>}"
threads="${THREADS:-16}"

deacon index build -k 31 -w 15 --threads "${threads}" --output "${output_index}" "${vpd_fasta}"
deacon index info "${output_index}"
