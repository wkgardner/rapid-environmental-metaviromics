# PLAN: v1 Rapid Wastewater Metagenomics pipeline

Source of truth: the approved spec, snapshot in `docs/v1-spec.md` (Claude Doc rev 377,
https://claude.ai/artifact/G491Kzg6ingETFz7t4517w). Vocabulary: `GLOSSARY.md`.
Settings rule: `docs/adr/0001-determined-defaults-exposed-params.md`.

## Ground rules (user, 2026-10-29)

- TDD, red then green, one vertical slice at a time.
- Small logical commits directly on `main`, detailed messages, no co-author or session lines.
- Stage files by explicit path only (no `.gitignore` yet; the user owns it).
- Local Nextflow execution only. No OSPool profile work for now.
- References and data live locally under project home (CHTC `/staging` is broken).
- Python via pixi only. No writes to `/tmp`; tools get `TMPDIR` under the work dir.

## Test seams (agreed)

| Seam | Tool | Catches | Misses |
| --- | --- | --- | --- |
| A. `bin/` scripts | pytest | table logic, LabKey TSV contract, spec edge cases | tool flags, channel wiring |
| B. each process | nf-test, real tool, tiny data | wrong flags, missing outputs | how processes chain |
| C. whole pipeline | nf-test on seeded synthetic data | wiring, skips, published layout | real-data scale, memory |

Fixtures for seam A come from real tool output already on disk
(`/path/to/rww-refbuild-scratch/demo_strain/`, `esv_masked_check/`), never invented formats.

## Layout

```
main.nf                      entry: parse params, call workflow, publish outputs (output {})
nextflow.config              params defaults, profiles (local is the only one exercised)
nextflow_schema.json         nf-schema param validation
assets/schema_input.json     samplesheet schema (NVD format)
workflows/rapid_ww.nf        wires subworkflows
subworkflows/local/          input_check, preprocessing, kmer_calling, mapping_calling, reporting
modules/local/<tool>/        one process per step
bin/                         verb-named Python (polars) scripts with all the logic
tests/python/                pytest for bin/ (seam A) + fixtures/
tests/                       nf-test for modules (seam B) and pipeline (seam C)
tests/data/                  output of `pixi run make-test-data` (seeded, tiny)
pixi.toml                    one env: nextflow, nf-test, every tool, python libs
```

## Tool pins (local pixi env)

nextflow 26.10.4 · nf-test 0.9.5 · deacon 0.18.0 (spec wants 0.19.1, not on bioconda yet) ·
bbmap >=39.19,<40 (NVD's pin) · sylph 1.0.0 · esviritu 1.3.5 · multiqc 1.35 · sra-tools 3.4.1 ·
python 3.12 + polars + biopython + pyyaml + pytest.

## Build order (each step = one or more red/green commits)

1. Repo bootstrap: docs snapshot, pixi env, PLAN/PROGRESS.
2. Test-data generator (`bin/make_test_data.py`, YAML spec, seeded): paired FASTQ + 7-record VPD subset,
   cases both-callers, k-mer-only, mapping-only, neither, empty sample.
3. Samplesheet: `bin/check_samplesheet.py` (NVD columns) + nf-schema + input channel; SRA fetch lifted from NVD.
4. Host depletion: Deacon `filter --deplete` process + summary JSON parsing.
5. Read QC: BBDuk trim, Clumpify dedup, BBDuk filter (NVD settings) + log parsing.
6. `preprocessing_read_counts` table (incl. `skipped` columns).
7. Genome k-mer-based calling: sylph sketch + profile processes, `bin/build_vpd_lookup.py`,
   `bin/build_genome_kmer_based_calls.py`.
8. Genome mapping-based calling: EsViritu process + wrapper status (`no_detections` / `failed`),
   `bin/build_genome_mapping_based_calls.py` (call_level, genus-name normalisation).
9. `sample_status`, `run_metadata.yaml`, MultiQC, published layout (01..05 folders).
10. End-to-end nf-test (seam C) on synthetic data.
11. Real-data run locally with user-supplied references; record timings.
12. Later: reference-build workflow, containers / mono image, OSPool profile.

## Open items

- Deacon 0.18.0 vs 0.19.1: confirm the user's union index works with 0.18.0.
- Resolved 2026-10-29: BBDuk counts `minlength` drops inside "Low quality discards" (254 = 200 low-quality + 54 short on the test data), so `qc_low_quality_reads_removed` covers quality and length.
