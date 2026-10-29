# Rapid Detection Wastewater Metagenomics v1 Spec

2026-10-28

v1 is an nf-core-style Nextflow pipeline for weekly wastewater viral surveillance with two callers: genome k-mer-based calling (sylph) reports first, and genome mapping-based calling (EsViritu) follows. Both search the same vector-masked VPD, and every output is an append-only TSV that a LabKey file watcher can import.

## Subway Diagram *W**H**E**E**E**E**E**E**!**!**!**!**!**!*

_[diagram: subway map, see the Doc]_

The input runs host depletion and read QC, then feeds sylph and EsViritu in parallel. Each caller writes its own table; comparing the two happens in SQL in the LIMS backend, not in the pipeline. The sylph path never waits for EsViritu.

## Inputs

v1 takes NVD's samplesheet unchanged and accepts SRA accessions, local FASTQ and ONT reads; every platform follows every step, handled the way NVD handles it.

- **Samplesheet:** `sample_id,srr,platform,fastq1,fastq2`, plus optional `fastq1_glob,fastq2_glob` for lane-split or barcode folders. Validated with nf-schema. No site or date metadata in v1.
- **SRA fetch:** lifted from NVD v3.5.
- **Read layouts:** Illumina paired, interleaved or single-end, and ONT single reads, as NVD supports them.
- **ONT is not a side branch:** ONT samples run host depletion, read QC, genome k-mer-based calling and genome mapping-based calling like any other sample. Only tool settings differ by platform: single-end Deacon, no adapter trimming and minimum average quality 12 in read QC, EsViritu `-mmP lr:hq`, sylph `-r`.

## Preprocessing

One unioned Deacon pass removes host and rRNA reads before read QC.

| Step | Tool | Input | Default settings | Output |
| --- | --- | --- | --- | --- |
| Host depletion | Deacon 0.19.1, `filter --deplete` | input reads | Unioned k31/w15 index: panhuman-1 (HPRC Year 1 ∪ CHM13v2.0 ∪ GRCh38.p14, minus FDA-ARGOS bacteria and RefSeq viruses) ∪ SILVA SSU 144 + LSU 138.2 NR99 (U→T). `abs 2`, `rel 0.01`, `complexity_threshold: null` | depleted reads |
| Read QC | BBTools 39.x: bbduk + clumpify (NVD's `modules/bbmap.nf`) | depleted reads | Three processes, NVD's settings. 1. Adapter trim (Illumina only): `bbduk.sh ref=adapters ktrim=r k=23 mink=11 hdist=1`, plus `tpe tbo removeifeitherbad=t` for pairs. 2. Dedup: `clumpify.sh dedupe=t subs=0 allowns=f optical=f reorder=f` (exact duplicates only). 3. Filter: `bbduk.sh minavgquality=20` (ONT 12) `minlength=50 entropy=0.5 entropywindow=50 entropyk=5`. | preprocessed reads → sylph + EsViritu |

- No separate telomere pass: telomeric reads pass depletion, but 2,000 simulated reads gave 0 sylph and 0 EsViritu calls.
- No tobamovirus depletion: EsViritu reports PMMoV and ToBRFV as ordinary calls.
- Read QC runs only on the main path; all read QC goes into one MultiQC report per run.

## Genome k-mer-based calling

A genome k-mer-based call is any sylph 1.0.0 profile row against the vector-masked VPD, using very-low-false-positive-rate settings similar to silly-fast-bfx and locked for v1.

```
sylph:
  database:      {kmer_length: 31, compression: 50}   # min_spacing: sylph default
  sample_sketch: {kmer_length: 31, compression: 50}
  profile:       {minimum_ani: 90, min_count_correct: 3, min_reference_kmers: 1,
                  min_matching_kmers: 0, redundancy_threshold: 99}   # -R, exposed
```

- The database is one FASTA per VPD Assembly, split at N stretches (sylph encodes N as A); segmented viruses count as one genome.
- ANI and abundance are columns, not gates. One threshold for all species.
- Its job is specificity: no false positives. Accepted in v1: a lentiviral vector calls HIV-1.
- Genome → VPD Assembly → species comes from the VPD metadata, joined on the sanitised `Genome_file` (123 Assembly names contain `/`).

## Genome mapping-based calling

A genome mapping-based call is an EsViritu 1.3.5 assembly with read_count ≥ 1, run on every sample against the same vector-masked VPD.

- Read from `detected_virus.assembly_summary.tsv` plus `info.tsv`. Breadth, depth, identity and RPKMF are columns, not gates. EsViritu `--dedup` is off (dedup already ran in read QC).
- **Genus-level calls:** below 0.90 identity EsViritu writes `s__unclassified_<genus>`. That row is `call_level = genus`, and everything below genus is null.
- **Wrapper:** EsViritu's exit codes are unreliable. "No reads aligned" (exit 1) → `no_detections`, a real negative; a crash or exit 0 with no tables → `failed`.
- Version 1.3.5 contains the fix for issue #51 (one SARS-CoV-2-rich sample needed 156 GiB on 1.3.4).

**Abundance normalisation (in v1):** every genome mapping-based call gets `normalized_read_ratio`, its `read_count` divided by the sample's total reads on a chosen genus. The genus is a param, `normalization_genus`, default `Tobamovirus`, matched by name on VPD's own `genus` metadata (`vpd_genus`).

- The denominator is the sum of `read_count` over that sample's genome mapping-based calls with `vpd_genus` = the chosen genus, at any `call_level` (species and genus-only rows both count). Dividing by it cancels sequencing depth and the compositional swings that rRNA and other abundant taxa cause.
- Why the genus, not PMMoV alone: a collaborator's July 2026 medRxiv preprint (not yet peer reviewed) found Tobamovirus-genus reads the best of seven denominators, improving median site-level concordance in 18 of 19 pathogen × comparison-source pairs (SARS-CoV-2 vs WW-PCR: Spearman R 0.38 → 0.72). ToBRFV was ~10× more abundant than PMMoV and ~72 % of genus reads.
- If a sample has no genome mapping-based call in the chosen genus, `normalized_read_ratio` is null for that sample.
- Normalisation applies to genome mapping-based calls only; genome k-mer-based calls are not normalised.

## Outputs

Every table is a TSV, one file per sample, written once when its step finishes and never rewritten, because LabKey's file drop can't update rows.

- **File name:** `<table>/<run_id>_<sample_id>.<table>.tsv`; an empty cell = null. Every row carries `run_id` (param; default `rww_<YYYYMMDD>_<Nextflow runName>`) and `sample_id`.
- **Tool-neutral names:** the tool and version behind each step live in `run_metadata.yaml`, so a tool can be swapped without changing the tables.

```
01_preprocessing/                  preprocessing_read_counts, Deacon summaries, BBDuk and Clumpify logs
02_genome_kmer_based_calling/      genome_kmer_based_calls
03_genome_mapping_based_calling/   genome_mapping_based_calls (+ EsViritu per-sample reports)
04_qc_report/                      <run_id>.multiqc_report.html
05_run_summary/                    sample_status, <run_id>.run_metadata.yaml
```

| Table | Grain | Published when | Columns (after `run_id`, `sample_id`) |
| --- | --- | --- | --- |
| `preprocessing_read_counts` | sample | read QC finishes | `raw_reads`, `host_reads_removed`, `reads_post_depletion`, `qc_adapter_trimmed_reads`, `qc_duplicate_reads_removed`, `qc_low_quality_reads_removed`, `qc_low_complexity_reads_removed`, `reads_post_qc` |
| `genome_kmer_based_calls` | sample × VPD Assembly | genome k-mer-based calling finishes | `vpd_assembly`, `vpd_species`, `vpd_subspecies`, `vpd_genus`, `vpd_family`, `vpd_taxids`, `taxonomic_abundance`, `sequence_abundance`, `adjusted_ani`, `naive_ani`, `effective_coverage`, `containment_matching_kmers`, `containment_total_kmers`, `kmers_reassigned` |
| `genome_mapping_based_calls` | sample × EsViritu Assembly | genome mapping-based calling finishes | `vpd_assembly`, `call_level`, `vpd_species`, `vpd_subspecies`, `vpd_genus`, `vpd_family`, `vpd_taxids`, `subspecies_label`, `read_count`, `covered_bases`, `assembly_length`, `breadth_fraction`, `mean_coverage`, `avg_read_identity`, `consensus_ref_identity`, `rpkmf`, `normalized_read_ratio` |
| `sample_status` | sample × step | each step ends (one file per step) | `step` (`preprocessing`, `genome_kmer_based_calling`, `genome_mapping_based_calling`), `status` (`ok`, `no_reads_left`, `skipped_no_reads`, `no_detections`, `failed`) |

- Read counts are reads (mates), not pairs. `host_reads_removed` is the one unioned pass (human + rRNA). `qc_adapter_trimmed_reads` counts reads BBDuk trimmed (kept, not removed). `qc_duplicate_reads_removed` is Clumpify's exact-duplicate count. `reads_post_qc` is EsViritu's RPKMF denominator.
- `vpd_taxids` is `;`-joined. `vpd_subspecies` is VPD's label for the reference (e.g. `SARS-CoV-2_Delta`): the nearest reference, not a variant call. `subspecies_label` is EsViritu's own adjusted label.
- `run_metadata.yaml`: run_id, pipeline version, start time, samplesheet, VPD/mask/taxdump versions, `steps:` (tool + version per step, each depletion pass with what it removes and its index), and every effective param.
- A throwaway prototype with fake data and a click-through viewer is in `.scratch/v1-pipeline-spec/prototypes/09-output-tables-PROTOTYPE/`.

### What the tables look like

Example rows from the prototype's fake run `rww_20261028_happy_turing` (real VPD names, fake counts). Each table shows its key columns; `run_id` is left off to save width.

**`preprocessing_read_counts`** (one row per sample)

| sample_id | raw_reads | host_reads_removed | reads_post_depletion | qc_duplicate_reads_removed | reads_post_qc |
| --- | --- | --- | --- | --- | --- |
| ww_madison_2026w43 | 4000000 | 1401204 | 2598796 | 290000 | 2221996 |
| ww_milwaukee_2026w43 | 6200000 | 2544310 | 3655690 | 410000 | 3132890 |
| negative_control_2026w43 | 12000 | 12000 | 0 |  | 0 |

**`genome_kmer_base``d_calls`** (one row per sylph hit; ww_madison has 11, four shown)

| sample_id | vpd_assembly | vpd_species | vpd_subspecies | adjusted_ani | effective_coverage | containment_matching_kmers | containment_total_kmers |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ww_madison_2026w43 | acc:MW423686.2 | Betacoronavirus pandemicum | SARS-CoV-2_Alpha | 99.61 | 3.9 | 180 | 183 |
| ww_madison_2026w43 | acc:BS001637.1 | Betacoronavirus pandemicum | SARS-CoV-2_Delta | 99.55 | 3.1 | 178 | 183 |
| ww_madison_2026w43 | acc:BS002303.1 | Betacoronavirus pandemicum | SARS-CoV-2_Omicron | 99.48 | 2.5 | 174 | 183 |
| ww_madison_2026w43 | set:EF564601 | Human mastadenovirus B | Human adenovirus 3+11p | 91.40 | 0.2 | 18 | 439 |

**`g``en``ome_mapping_based``_calls`** (ww_milwaukee; PMMoV is this sample's only Tobamovirus call, so the genus total is 70,114 reads and PMMoV's `normalized_read_ratio` is 1)

| sample_id | vpd_assembly | call_level | vpd_species | vpd_genus | read_count | breadth_fraction | rpkmf | normalized_read_ratio |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ww_milwaukee_2026w43 | GCF_000859645.1 | species | Tobamovirus capsici | Tobamovirus | 70114 | 1.000 | 3520.5 | 1 |
| ww_milwaukee_2026w43 | acc:MW423686.2 | species | Betacoronavirus pandemicum | Betacoronavirus | 1204 | 1.000 | 12.88 | 0.0172 |
| ww_milwaukee_2026w43 | set:84ac15bd… | species | Alphainfluenzavirus influenzae | Alphainfluenzavirus | 812 | 1.000 | 19.44 | 0.0116 |
| ww_milwaukee_2026w43 | GCF_000846685.1 | genus |  | Mastadenovirus | 388 | 1.000 | 3.62 | 0.00553 |
| ww_milwaukee_2026w43 | GCF_000857045.1 | species | Orthopoxvirus monkeypox | Orthopoxvirus | 2 | 0.002 | 0.0032 | 0.0000285 |

**`sample_status`** (one file per step)

| sample_id | step | status |
| --- | --- | --- |
| ww_lacrosse_2026w43 | preprocessing | ok |
| ww_lacrosse_2026w43 | genome_kmer_based_calling | ok |
| ww_lacrosse_2026w43 | genome_mapping_based_calling | no_detections |
| ww_greenbay_2026w43 | genome_mapping_based_calling | failed |
| negative_control_2026w43 | genome_kmer_based_calling | skipped_no_reads |

Blank cells are nulls: no read QC on an empty sample, no species below a genus call.

## LIMS upload

The pipeline only writes files; getting them into the LIMS happens outside it.

- **LabKey:** file watchers on a LabKey container pick up the tables, delivered over WebDAV.
- **S3:** the outputs can also be uploaded directly to an S3 bucket.

## Parameters and defaults

Every tool setting and threshold has a spec-chosen default and is exposed as a param in the YAML params file (nf-schema validated); nothing is hard-coded in modules or `bin/` (ADR 0001).

| Group | Defaults |
| --- | --- |
| Host depletion | `abs_threshold 2`, `rel_threshold 0.01`, `complexity_threshold null` |
| Read QC (BBTools) | `trim_adapters true` (Illumina only), `dedup true` (Clumpify, exact), `min_read_quality_illumina 20`, `min_read_quality_nanopore 12`, `min_read_length 50`, `filter_low_complexity_reads true`, `min_read_entropy 0.5` (NVD's names and values) |
| sylph | k31, `-c 50`, `minimum_ani 90`, `min_count_correct 3`, `min_reference_kmers 1`, `min_matching_kmers 0`, `redundancy_threshold 99` (locked) |
| Genome mapping-based calls | `mapping_based_call_min_read_count 1`, EsViritu `dedup false` |
| Normalisation | `normalization_genus Tobamovirus` (matched on VPD's genus name) |
| Run | `run_id` (default `rww_<YYYYMMDD>_<runName>`), `mono_container false` |

Changing a default means changing the spec. A per-run override goes in that run's params YAML, and the effective values land in `run_metadata.yaml`. All configs are YAML only.

## References, tool versions and containers

Both callers search VPD v3.2.4 with a vector-only mask: 13.7 kb (0.0058 %) N-masked, which removed all four simulated cloning-vector false calls with zero true-call loss.All references will be included as build processes so anyone can download the same references.

- **Mask sources:** UniVec_Core, FCS-adaptor, lab vectors and NAO contaminants, kept as a YAML manifest (accession, source, reason) versioned with VPD. No SILVA mask on VPD (Deacon handles rRNA).
- **Indexes:** union depletion index (~3.5 GB, ~4.5 GB RAM to filter), sylph VPD database (k31 c50), EsViritu DB (447 MB tarball, ~894 MB unpacked). A reference-build workflow, run once per VPD and mask version, produces and stages them.

| Tool | Version | Container (per-process mode) |
| --- | --- | --- |
| Deacon | 0.19.1 (30 % faster gzip) | `quay.io/biocontainers/deacon:0.19.1` once bioconda has it; else our own |
| BBTools (bbduk, clumpify) | 39.x (NVD allows ≥39.19, <40) | bioconda bbmap; processes lifted from NVD `modules/bbmap.nf` |
| sylph / sylph-tax | 1.0.0 / 1.9.2 | nf-core modules, version bump |
| EsViritu | 1.3.5 | biocontainers, local module behind a wrapper |
| `bin/` Python | polars, biopython, pyyaml | one small image built from a pixi env |
| All of the above | same pins | **mono image** (one pixi env with every tool), used when `mono_container: true` |

## Engineering conventions

Build a standard nf-core-style pipeline first, then make it run on OSPool through a profile.

- **Nextflow:** newest release, strict syntax, workflow outputs (`output {}`). No typed processes in v1: ordinary untyped processes and plain meta maps (static typing is a preview nf-core hasn't adopted).
- **Layout:** nf-core-style folders, every step its own process.
- **Skipping steps:** every process has a `--skip_<process>` param (for example `--skip_host_depletion`, `--skip_dedup`), default false. A skipped step passes its input straight through, so the rest of the pipeline still runs. Its columns in the output tables read `skipped` instead of a count (for example `qc_duplicate_reads_removed = skipped`), and `sample_status` records the step with status `skipped`.
- **Containers:** per-process containers by default. The flag `mono_container: true` runs every process in one image instead (one pixi env with every pinned tool), set through a config profile so modules stay unchanged.
- **AWS cloud** execution support.
- **Minimal Groovy:** `.nf` files only wire channels. All logic is in verb-named Python scripts in `bin/`, using polars, not hand-rolled CSV. Tests: pytest + nf-test. Environments: pixi.
- **Module reuse order:** nf-core, then silly-fast-bfx, then NVD (lift specific modules only).
- **Read filtering:** Deacon for host depletion; BBTools (bbduk, clumpify) for read QC, lifted from NVD so dedup matches NVD exactly.
- **OSPool:** nf-ospool v0.1.1 ignores `container_image`, so the `ospool` profile sets it through `clusterOptions`. Large references are plain path params under a versioned `/staging/groups/<group>/refs/`.
- **Test data:** a seeded Biopython generator driven by a YAML spec (`pixi run make-test-data`) writes paired FASTQ and a 7-record VPD subset (<1 MB, no network), covering both-callers, k-mer-only, mapping-only and neither cases.

## Watchlist

The watchlist is a YAML file of categories; v1 ships only this skeleton, and a colleague builds out the categories.

```
watchlist_version: "0.0.0-skeleton"
categories:
  - name: example_category   # unique, snake_case
    priority: 1              # integer; higher = more severe (surveillance level, not biosafety level)
    include: [0000000]       # NCBI TaxIDs; descendants included
    exclude: []              # NCBI TaxIDs removed from include
    rg_human: null           # optional WHO/PHAC risk group 1-4
    notes: ""
    sources: []              # e.g. Tisza ontology, Jaffe risk groups
```

- Resolved once per VPD version, against a pinned NCBI taxdump, into a VPD Assembly → category, priority table.
- Human pathogens only in v1. Variola gets its own high-priority category.

## Out of scope and future release ideas

v1 leaves comparing the two callers' calls, rRNA bait screening and LIMS upload outside the pipeline.

**Out of scope for v1**

- Any strain-level comparison or corroborate/override of genome k-mer-based and genome mapping-based calls, in the pipeline or in post-processing.
- Comparing the two callers (for example species confirmation): done in SQL in the LIMS backend over `genome_kmer_based_calls` and `genome_mapping_based_calls`.
- rRNA extraction and bait screening: no rRNA branch in v1; silly-fast-bfx stays a separate tool.
- Uploading to LabKey or any LIMS: file watchers over WebDAV, or a direct S3 upload (see LIMS upload).
- A SQLite-watch trigger for launching runs: v1 starts from a samplesheet. Add in sqlite db watcher as our sample manifest from sluice is built for auto kickoffs.

**Future release ideas**

- Risk-group categorization of calls: join the watchlist's Assembly → category, priority table onto the call tables (shown dashed on the subway map).
- An agent process that builds weather-style reports: readable weekly summaries and longitudinal trends from the call tables.
- Typed Nextflow processes and records, once nf-core supports them.
- Watchlist-specific sylph thresholds.
- Mask lentiviral-vector backbones in VPD, after measuring the HIV-1 sensitivity cost.

## Open implementation tasks

The two biggest risks are EsViritu's memory on real samples and OSPool queue latency against the same-day target.

- [ ] Resource profiling for an average SB weekly plate.
- [ ] Get the Deacon 0.19.1 container (bioconda bump, or our own image).
- [ ] Build the mono image (one pixi env, every pinned tool) and the config profile behind `mono_container`.
- [ ] Check ONT end to end: ONT reads in the test generator, Clumpify and BBDuk on long reads (NVD runs both on every platform), and sylph/EsViritu behaviour on long reads (thresholds were measured on Illumina).
- [ ] Pin the exact BBTools version (NVD allows 39.19 to <40) and check which BBDuk log line counts `minlength` drops, then map BBDuk/Clumpify logs onto the `qc_*` columns and into MultiQC.
- [ ] Shape the reference-build workflow and where its outputs are staged. Make it reusable anf fully transparent for other people who want to use the same references. Like SB did!
- [ ] Set up the repo skeleton and CI (GitHub Actions, nf-test and pytest layout).
- [ ] Build the watchlist categories, converting the Tisza and Jaffe lists and checking the possible Tai Forest/Zaire swap in the Jaffe table.
- [ ] Not yet specified: typing below species (flu serotypes, mpox clades), the history store for trends, the weekly "weather report" agent, and an optional core_nt BLAST check for priority calls.
