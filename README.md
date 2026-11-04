# Rapid Environmental MetaViromics (REMV)

REMV is a Nextflow pipeline for weekly viral surveillance from wastewater sequencing. It cleans the reads, then calls viruses two ways against the same reference database, EsViritu's Virus Pathogen Database (VPD):

- **Genome k-mer-based calling** (sylph) is fast and tuned for a very low false-positive rate.
- **Genome mapping-based calling** (EsViritu) maps every read, so it is slower but finds more low-abundance viruses.

Every result is a per-sample TSV that is written once and never rewritten, so a LIMS file watcher can import the tables as they appear.

![Pipeline diagram](docs/metro_map.svg)

> **Status: 1.0.0dev.** The pipeline runs locally through a pixi environment on `linux-64`. Per-process containers and an OSPool profile are planned but not built yet.

## Quick start

You need [pixi](https://pixi.sh) and the references listed below.

```bash
git clone https://github.com/wkgardner/rapid-environmental-metagenomics.git
cd rapid-environmental-metagenomics
pixi install

cp conf/params_template.yaml my_run.yaml   # then fill in your paths
pixi run nextflow run . -params-file my_run.yaml
```

On a shared machine, add `-c conf/shared_host.config` to cap the run at 16 CPUs and 160 GB, and use `-work-dir` to keep the work directory somewhere with enough space.

## Samplesheet

The pipeline uses [NVD](https://github.com/dholab/nvd)'s samplesheet format, with one row per sample:

```csv
sample_id,srr,platform,fastq1,fastq2,fastq1_glob,fastq2_glob
site_a_week44,,illumina,/data/site_a_R1.fastq.gz,/data/site_a_R2.fastq.gz,,
site_b_week44,,illumina,,,/data/site_b_L*_R1_001.fastq.gz,/data/site_b_L*_R2_001.fastq.gz
site_c_week44,,ont,/data/site_c.fastq.gz,,,
public_sample,SRR33296246,illumina,,,,
```

- `platform` is required: `illumina` or `ont`.
- FASTQ paths and glob patterns must be absolute.
- Use the glob columns for samples split across lanes.

To build a samplesheet from a folder of bcl2fastq or BCL Convert files:

```bash
pixi run python bin/write_samplesheet_from_fastq_dir.py --fastq-dir /data/run --out samplesheet.csv
```

## References

| Param | What it is | How to get it |
| --- | --- | --- |
| `reference_filter_index` | Deacon index of the reads to remove (default mode) or keep | `scripts/build_host_depletion_index.sh` builds panhuman-1 + SILVA SSU/LSU, after you download the inputs listed at the top of the script |
| `esviritu_database` | Folder holding `virus_pathogen_database.{fna,fna.fai,mmi,all_metadata.tsv}` | The [EsViritu database](https://doi.org/10.5281/zenodo.17716199) |
| `vpd_lookup` | One row per VPD assembly, shared by both call tables | `bin/build_vpd_lookup.py` |
| `sylph_database` | VPD sketched at k31, c50 | `bin/write_sylph_genome_fastas.py`, then `sylph sketch` |

To build the last two from a VPD folder (`$VPD`):

```bash
pixi run python bin/build_vpd_lookup.py \
    --vpd-metadata $VPD/virus_pathogen_database.all_metadata.tsv --out vpd_lookup.tsv
pixi run python bin/write_sylph_genome_fastas.py \
    --vpd-fasta $VPD/virus_pathogen_database.fna \
    --vpd-metadata $VPD/virus_pathogen_database.all_metadata.tsv \
    --vpd-lookup vpd_lookup.tsv --min-fragment-length 31 \
    --out-dir sylph_genomes --genome-list sylph_genome_list.txt
pixi run sylph sketch --gl sylph_genome_list.txt -k 31 -c 50 -o vpd
```

> **Not reproducible from this repo yet:** our production runs use VPD 3.2.4 with a vector-only mask and a reference filter index that also includes Rfam 15.1 rRNA families. The build steps for the mask and for the Rfam addition are not scripted here, so the commands above give you an unmasked VPD and a filter index without Rfam.

If your filter index differs from the default, set `reference_filter_index_description` so that `run_metadata.yaml` describes it correctly. The template also sets `vpd_version` and `vpd_mask_version` for the same reason.

## Outputs

```
<outdir>/
├── 01_preprocessing/
│   ├── preprocessing_read_counts/
│   └── logs/
├── 02_genome_kmer_based_calling/
│   └── genome_kmer_based_calls/
├── 03_genome_mapping_based_calling/
│   ├── genome_mapping_based_calls/
│   └── esviritu_reports/<sample_id>/esviritu_output/   (and <sample_id>_temp/ with esviritu_keep)
├── 04_qc_report/                      <run_id>.multiqc_report.html
└── 05_run_summary/
    ├── sample_status/
    └── <run_id>.run_metadata.yaml
```

Tables are named `<run_id>_<sample_id>.<table>.tsv`, and every row carries `run_id` and `sample_id`. An empty cell means null.

| Table | One row per | Key columns |
| --- | --- | --- |
| `preprocessing_read_counts` | sample | `raw_reads`, `reads_removed_by_reference_filter`, `reads_post_reference_filter`, `qc_duplicate_reads_removed`, `reads_post_qc` |
| `genome_kmer_based_calls` | sample × VPD assembly | `vpd_assembly`, `vpd_species`, `vpd_genus`, `adjusted_ani`, `effective_coverage`, `taxonomic_abundance` |
| `genome_mapping_based_calls` | sample × VPD assembly | `vpd_assembly`, `call_level`, `vpd_species`, `vpd_genus`, `read_count`, `breadth_fraction`, `rpkmf`, `normalized_read_ratio` |
| `sample_status` | sample × step | `step`, `status` (`ok`, `no_reads_left`, `skipped_no_reads`, `no_detections`, `failed`, `skipped`) |

`run_metadata.yaml` records the tool and version behind each step, the reference versions and every effective param. Table and column names stay tool-neutral, so a tool can change without changing the tables.

### Reading the calls

- **The two callers disagree at the edges, by design.** EsViritu reports more species, many of them on a handful of reads. sylph can miss viruses present at low levels in large genomes.
- **The pipeline does not compare the callers.** Compare `genome_kmer_based_calls` with `genome_mapping_based_calls` downstream, for example in SQL in your LIMS. See [GLOSSARY.md](GLOSSARY.md) for the terms.
- **`vpd_subspecies` is not a strain call.** It names the nearest reference genome, not a variant detected in the sample.
- **`normalized_read_ratio` is relative to one genus.** It divides a call's `read_count` by the sample's total reads in `normalization_genus` (default `Tobamovirus`). It is null when the sample has no call in that genus, and it applies to mapping-based calls only.

## Configuration

Every tool setting has a default in [nextflow.config](nextflow.config) and can be overridden in your params YAML ([why](docs/adr/0001-determined-defaults-exposed-params.md)). The most common switches are:

| Param | Default | Effect |
| --- | --- | --- |
| `reference_filter_mode` | `deplete` | `deplete` removes reads matching the index; `extract` keeps only those reads |
| `skip_reference_filter` | `false` | Skips Deacon; the post-filter read count is left null |
| `trim_adapters` | `true` | BBDuk adapter trimming (Illumina only) |
| `dedup` | `true` | Clumpify exact-duplicate removal |
| `skip_read_filtering` | `false` | Skips the BBDuk quality, length and low-complexity filter |
| `skip_genome_kmer_based_calling` | `false` | Skips sylph; `sylph_database` is then not needed (`vpd_lookup` still is) |
| `skip_genome_mapping_based_calling` | `false` | Skips EsViritu; `esviritu_database` is then not needed (`vpd_lookup` still is) |
| `esviritu_keep` | `true` | Publishes EsViritu's temp folder, including its BAM files, beside its raw files. These add up across a run; set `false` to save space |
| `normalization_genus` | `Tobamovirus` | Denominator for `normalized_read_ratio` |

The sylph settings (k31, c50, minimum ANI 90) are the v1 defaults. The sylph database must be built at the same k-mer length and compression as the sample sketches.

## Tests

```bash
pixi run test-python      # bin/ scripts (pytest)
pixi run test-nf          # modules and the whole pipeline on synthetic data (nf-test)
pixi run make-test-data   # regenerate the synthetic data
```

## Tools

[Deacon](https://github.com/bede/deacon), [BBTools](https://bbmap.org/), [sylph](https://github.com/bluenote-1577/sylph), [EsViritu](https://github.com/cmmr/EsViritu), [MultiQC](https://multiqc.info/), [sra-tools](https://github.com/ncbi/sra-tools), [sracha](https://github.com/rnabioco/sracha-rs) and [Nextflow](https://www.nextflow.io/). Exact versions are pinned in [pixi.toml](pixi.toml). Please cite these tools if you use the pipeline.

## License

MIT
