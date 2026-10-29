# PROGRESS

See `PLAN.md` for the build order. Update after every step.

## Status

- [x] 1. Repo bootstrap (spec snapshot, pixi env, plan)
- [x] 2. Test-data generator + committed synthetic dataset (`pixi run make-test-data`)
- [x] 3. Samplesheet: NVD read-input code lifted (`bin/resolve_read_inputs.py`), RESOLVE_READ_INPUTS
- [x] 4. Host depletion: HOST_DEPLETION (local) and HOST_DEPLETION_SRA (sracha | deacon)
- [x] 5. Read QC: BBDUK_TRIM_ADAPTERS, CLUMPIFY_DEDUP, BBDUK_FILTER_READS
- [x] 6. preprocessing_read_counts builder
- [x] 7. Genome k-mer-based calling: VPD lookup, sylph genomes, sketch, profile, call table
- [x] 8. Genome mapping-based calling: ESVIRITU wrapper, call table with genus normalisation
- [x] 9. sample_status, run_metadata.yaml, MultiQC (04_qc_report), published layout
- [x] 10. Workflow wiring (`workflows/rapid_ww.nf`, `main.nf`) + end-to-end nf-test
- [ ] 11. Real-data local run: smoke run on SRR33296246 with skip_host_depletion (no union index yet)
- [ ] 12. nf-schema validation (`nextflow_schema.json`, samplesheet schema) — not started
- [ ] 13. Reference-build workflow, containers / mono image, OSPool profile — later

## How to test

- `pixi run test-python` (seam A, pytest: 50 tests)
- `pixi run nf-test test modules/local` (seam B: 16 tests; HOST_DEPLETION_SRA is tagged `network`)
- `pixi run nf-test test tests/main.nf.test` (seam C: whole pipeline on the synthetic data, ~11 s)

## Real references (local, /staging is down)

- sylph DB + VPD lookup: `.scratch/refs/vpd_3.2.4_vector_only/` (built 2026-10-29 with the
  pipeline's scripts: 19,925 assemblies, 38 MB k31/c50 .syldb, 6 s)
- EsViritu DB: `/path/to/rww-refbuild-scratch/esv_masked_check/db/v3.2.4_vector_only`
- Deacon union index: the user is building it.

## Log

- 2026-10-29: spec approved (Doc rev 377); plan written.
- 2026-10-29: steps 1–10 done. Findings: sracha 0.8 prints 15 info columns (NVD's 8-column check
  would fail; we read by column name). BBDuk counts minlength drops in "Low quality discards".
  The generator now mutates a diverged source once (independent per-read errors let sylph
  recover reference k-mers). EsViritu relabels genus-level calls only in tax_profile.tsv, so the
  call-table builder applies EsViritu's rule. MultiQC 1.35 needs a search-pattern tweak for
  BBTools 39.x BBDuk logs. Spec Doc updated to rev 378 with these.

## Overnight plan (2026-10-29, user away)

1. Real run on two pilot samples: `runs/rww_20261029_real/` (params.yaml, run.out, trace.txt).
   First pilot sample done; the second's EsViritu running. Fix anything that breaks, then report timings/calls.
2. batch2: rsync of the selected sites into `real-samples/batch2-subset/` (8 lanes/sample, ~1 TB+).
   When no rsync into it remains: `bin/write_samplesheet_from_fastq_dir.py --fastq-dir
   real-samples/batch2-subset --out real-samples/batch2-subset/samplesheet.csv`, then launch
   `runs/rww_<date>_batch2/` with the same references and `-c conf/shared_host.config`
   (16 CPUs / 160 GB cap for the shared AP).
3. Tower: the user's `tower` block is in nextflow.config, uncommitted (do not commit it until the
   user confirms the token is not in the file).

## Next action

Read the smoke-run timings; run on real samples once the Deacon union index is ready.
