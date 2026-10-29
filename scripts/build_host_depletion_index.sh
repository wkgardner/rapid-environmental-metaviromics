#!/usr/bin/env bash
# Build the spec's unioned host-depletion index: panhuman-1 + SILVA SSU 144 + LSU 138.2 NR99,
# all k31/w15. SILVA is RNA, so U is converted to T first (spec). Run from the repo root:
#   pixi run bash scripts/build_host_depletion_index.sh ref
# Inputs under <ref_dir>: silva/SILVA_144_SSURef_NR99_tax_silva_trunc.fasta.gz,
# silva/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz, deacon-indexes/panhuman-1.k31w15.idx.
# Outputs in <ref_dir>/deacon-indexes/: silva_ssu144.k31w15.idx, silva_lsu138.2.k31w15.idx,
# silva_ssu144_lsu138.2.k31w15.idx and panhuman-1_silva_ssu144_lsu138.2.k31w15.idx.
set -euo pipefail

ref_dir="${1:?usage: build_host_depletion_index.sh <ref_dir>}"
threads="${THREADS:-16}"
deacon_dir="${ref_dir}/deacon-indexes"
export TMPDIR="${TMPDIR:-${ref_dir}/tmp}"
mkdir -p "${TMPDIR}"

build_rrna_index() {
    local silva_fasta="$1" index_out="$2"
    # U -> T, then index at the pipeline's k31/w15.
    seqkit seq --rna2dna --threads "${threads}" "${silva_fasta}" \
    | deacon index build -k 31 -w 15 --threads "${threads}" --output "${index_out}" -
}

build_rrna_index "${ref_dir}/silva/SILVA_144_SSURef_NR99_tax_silva_trunc.fasta.gz" "${deacon_dir}/silva_ssu144.k31w15.idx"
build_rrna_index "${ref_dir}/silva/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz" "${deacon_dir}/silva_lsu138.2.k31w15.idx"

deacon index union --output "${deacon_dir}/silva_ssu144_lsu138.2.k31w15.idx" \
    "${deacon_dir}/silva_ssu144.k31w15.idx" "${deacon_dir}/silva_lsu138.2.k31w15.idx"

deacon index info "${deacon_dir}/panhuman-1.k31w15.idx"
deacon index union --output "${deacon_dir}/panhuman-1_silva_ssu144_lsu138.2.k31w15.idx" \
    "${deacon_dir}/panhuman-1.k31w15.idx" "${deacon_dir}/silva_ssu144_lsu138.2.k31w15.idx"

for index in silva_ssu144 silva_lsu138.2 silva_ssu144_lsu138.2 panhuman-1_silva_ssu144_lsu138.2; do
    echo "== ${index}"
    deacon index info "${deacon_dir}/${index}.k31w15.idx"
done
