# PROGRESS

See `PLAN.md` for the build order. Update after every step.

## Status

- [x] 1. Repo bootstrap (spec snapshot, pixi env, plan)
- [x] 2. Test-data generator + committed synthetic dataset (`pixi run make-test-data`)
- [x] 3. Samplesheet: NVD read-input code lifted (`bin/resolve_read_inputs.py`), RESOLVE_READ_INPUTS
- [x] 4. Host depletion: HOST_DEPLETION (local) and HOST_DEPLETION_SRA (sracha | deacon)
- [x] 5. Read QC: BBDUK_TRIM_ADAPTERS, CLUMPIFY_DEDUP, BBDUK_FILTER_READS
- [x] 6. preprocessing_read_counts builder (`bin/build_preprocessing_read_counts.py`)
- [ ] 7. Genome k-mer-based calling: VPD lookup, sylph genomes, SYLPH_SKETCH_READS, SYLPH_PROFILE done;
        `bin/build_genome_kmer_based_calls.py` next
- [ ] 8. Genome mapping-based calling (EsViritu)
- [ ] 9. sample_status, run_metadata, MultiQC, publish layout
- [ ] 10. Workflow wiring + end-to-end nf-test
- [ ] 11. Real-data local run

## How to test

- `pixi run test-python` (seam A, pytest)
- `pixi run nf-test test modules/local` (seam B; the HOST_DEPLETION_SRA test is tagged `network`)

## Log

- 2026-10-29: spec approved (Doc rev 377); plan written.
- 2026-10-29: steps 1–6 done. Findings: sracha 0.8 prints 15 info columns (NVD's 8-column check
  would fail; we read by column name). BBDuk counts minlength drops in "Low quality discards".
  The generator now mutates a diverged source once (independent per-read errors let sylph
  recover reference k-mers and call a 12 %-diverged genome at 96 % ANI).

## Next action

TDD `bin/build_genome_kmer_based_calls.py` with a real sylph profile fixture, then EsViritu.
