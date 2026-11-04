// Rapid wastewater metagenomics v1: wiring only. All decisions about counts, calls and statuses
// live in the bin/ Python scripts the processes run.

include { RESOLVE_READ_INPUTS              } from '../modules/local/resolve_read_inputs'
include { REFERENCE_FILTER                 } from '../modules/local/reference_filter'
include { REFERENCE_FILTER_SRA             } from '../modules/local/reference_filter_sra'
include { PREPARE_READS                    } from '../modules/local/prepare_reads'
include { FETCH_SRA_READS                  } from '../modules/local/fetch_sra_reads'
include { BBDUK_TRIM_ADAPTERS              } from '../modules/local/bbduk_trim_adapters'
include { CLUMPIFY_DEDUP                   } from '../modules/local/clumpify_dedup'
include { BBDUK_FILTER_READS               } from '../modules/local/bbduk_filter_reads'
include { BUILD_PREPROCESSING_READ_COUNTS  } from '../modules/local/build_preprocessing_read_counts'
include { SYLPH_SKETCH_READS               } from '../modules/local/sylph_sketch_reads'
include { SYLPH_PROFILE                    } from '../modules/local/sylph_profile'
include { BUILD_GENOME_KMER_BASED_CALLS    } from '../modules/local/build_genome_kmer_based_calls'
include { ESVIRITU                         } from '../modules/local/esviritu'
include { BUILD_GENOME_MAPPING_BASED_CALLS } from '../modules/local/build_genome_mapping_based_calls'
include { WRITE_SAMPLE_STATUS              } from '../modules/local/write_sample_status'
include { MULTIQC                          } from '../modules/local/multiqc'

// How many read-QC logs a sample with reads left will produce, so its read-count row can be
// built as soon as its own logs exist instead of waiting for every sample.
def countReadQcLogs(meta) {
    def trims = params.trim_adapters && meta.platform == 'illumina' ? 1 : 0
    def dedups = params.dedup ? 1 : 0
    def filters = params.skip_read_filtering ? 0 : 1
    return trims + dedups + filters
}

workflow RAPID_WW {
    take:
    samplesheet     // path: NVD-format samplesheet
    run_id          // val: rww_<YYYYMMDD>_<runName> unless set

    main:
    // ---- Inputs: NVD's resolver writes one JSON record per sample ----
    RESOLVE_READ_INPUTS(samplesheet)
    ch_resolved = RESOLVE_READ_INPUTS.out.jsonl
        .splitText()
        .filter { line -> line.trim() }
        .map { line -> new groovy.json.JsonSlurper().parseText(line) }
        .branch { record ->
            sra: record.source == 'sra'
            local: true
        }
    ch_local_reads = ch_resolved.local.map { record ->
        def is_paired = record.read_structure == 'paired'
        def read_paths = is_paired ? record.r1 + record.r2 : record.reads
        def meta = [id: record.sample_id, platform: record.platform, read_mode: is_paired ? 'paired' : 'single',
                    r1_count: is_paired ? record.r1.size() : record.reads.size()]
        tuple(meta, read_paths.collect { read_path -> file(read_path) })
    }
    ch_sra_runs = ch_resolved.sra.map { record ->
        [id: record.sample_id, platform: record.platform, srr: record.srr.toUpperCase()]
    }

    // ---- Reference filter: one Deacon pass (deplete or extract); pairs come out interleaved ----
    // --skip_reference_filter skips it entirely: PREPARE_READS / FETCH_SRA_READS only merge lanes
    // and interleave pairs, and the read count after the filter is unknown (null), not guessed.
    if (params.skip_reference_filter) {
        ch_filtered = PREPARE_READS(ch_local_reads).reads
            .mix(FETCH_SRA_READS(ch_sra_runs).reads.map { meta, reads, read_mode -> tuple(meta + [read_mode: read_mode], reads) })
            .map { meta, reads -> tuple(meta + [reads_post_reference_filter: null], reads) }
        ch_reference_filter_summaries = channel.empty()
    } else {
        def filter_index = file(params.reference_filter_index, checkIfExists: true)
        ch_filtered = REFERENCE_FILTER(ch_local_reads, filter_index, params.reference_filter_mode).reads
            .mix(REFERENCE_FILTER_SRA(ch_sra_runs.map { meta -> tuple(meta, filter_index) }, params.reference_filter_mode).reads
                .map { meta, reads, read_count, read_mode -> tuple(meta + [read_mode: read_mode], reads, read_count) })
            .map { meta, reads, read_count -> tuple(meta + [reads_post_reference_filter: read_count as long], reads) }
        ch_reference_filter_summaries = REFERENCE_FILTER.out.summary.mix(REFERENCE_FILTER_SRA.out.summary)
    }
    ch_depleted = ch_filtered.branch { meta, reads ->
        has_reads: meta.reads_post_reference_filter == null || meta.reads_post_reference_filter > 0
        no_reads_left: true
    }

    // ---- Read QC: NVD's BBTools chain (trim -> dedup -> filter) ----
    ch_for_trimming = ch_depleted.has_reads.branch { meta, reads ->
        trim: params.trim_adapters && meta.platform == 'illumina'
        pass_through: true
    }
    BBDUK_TRIM_ADAPTERS(ch_for_trimming.trim)
    ch_trimmed = BBDUK_TRIM_ADAPTERS.out.reads.mix(ch_for_trimming.pass_through)

    ch_deduplicated = params.dedup ? CLUMPIFY_DEDUP(ch_trimmed).reads : ch_trimmed

    // reads_post_qc is the filter's output count, else the reference-filter count (null if unknown).
    ch_preprocessed = params.skip_read_filtering
        ? ch_deduplicated.map { meta, reads -> tuple(meta, reads, meta.reads_post_reference_filter) }
        : BBDUK_FILTER_READS(ch_deduplicated).reads
    ch_preprocessed_by_reads = ch_preprocessed
        .map { meta, reads, read_count -> tuple(meta + [reads_post_qc: read_count == null ? null : read_count as long], reads) }
        .branch { meta, reads ->
            callable: meta.reads_post_qc == null || meta.reads_post_qc > 0
            no_reads_left: true
        }

    // ---- preprocessing_read_counts: built per sample once its own logs exist ----
    ch_step_logs = ch_reference_filter_summaries
        .mix(BBDUK_TRIM_ADAPTERS.out.log)
        .mix(params.dedup ? CLUMPIFY_DEDUP.out.log : channel.empty())
        .mix(params.skip_read_filtering ? channel.empty() : BBDUK_FILTER_READS.out.log)
        .map { meta, log_file -> tuple(meta.id, log_file) }
    def reference_filter_logs = params.skip_reference_filter ? 0 : 1
    ch_expected_log_counts = ch_depleted.has_reads.map { meta, reads -> tuple(meta.id, meta, reference_filter_logs + countReadQcLogs(meta)) }
        .mix(ch_depleted.no_reads_left.map { meta, reads -> tuple(meta.id, meta, reference_filter_logs) })
    ch_sample_logs = ch_step_logs
        .combine(ch_expected_log_counts, by: 0)
        .map { sample_id, log_file, meta, expected_logs -> tuple(groupKey(sample_id, expected_logs), meta, log_file) }
        .groupTuple()
        .map { sample_key, metas, log_files -> tuple(metas[0], log_files) }
    BUILD_PREPROCESSING_READ_COUNTS(ch_sample_logs, run_id)

    // ---- Genome k-mer-based calling (sylph) ----
    def vpd_lookup = file(params.vpd_lookup, checkIfExists: true)
    ch_callable = ch_preprocessed_by_reads.callable
    ch_kmer_calls = params.skip_genome_kmer_based_calling ? channel.empty()
        : BUILD_GENOME_KMER_BASED_CALLS(
            SYLPH_PROFILE(SYLPH_SKETCH_READS(ch_callable).sketch, file(params.sylph_database, checkIfExists: true)).profile,
            vpd_lookup, run_id).table

    // ---- Genome mapping-based calling (EsViritu) ----
    ch_mapping_status = params.skip_genome_mapping_based_calling ? channel.empty()
        : ESVIRITU(ch_callable, file(params.esviritu_database, checkIfExists: true)).status
    ch_mapping_tables = params.skip_genome_mapping_based_calling ? channel.empty()
        : ESVIRITU.out.tables.map { meta, assembly_summary, info, status -> tuple(meta.id, [assembly_summary, info]) }
    ch_mapping_results = ch_mapping_status
        .map { meta, status -> tuple(meta.id, meta, status) }
        .join(ch_mapping_tables, by: 0, remainder: true)
        .map { sample_id, meta, status, tables -> tuple(meta, status, tables ?: []) }
    ch_mapping_calls = params.skip_genome_mapping_based_calling ? channel.empty()
        : BUILD_GENOME_MAPPING_BASED_CALLS(ch_mapping_results, vpd_lookup, run_id).table

    // ---- sample_status: one row per sample and step ----
    ch_samples_without_reads = ch_depleted.no_reads_left.map { meta, reads -> meta }
        .mix(ch_preprocessed_by_reads.no_reads_left.map { meta, reads -> meta })
    ch_status_rows = ch_preprocessed_by_reads.callable
        .map { meta, reads -> tuple(meta, 'preprocessing', meta.reads_post_qc == null ? [status: 'ok'] : [reads_left: meta.reads_post_qc], []) }
        .mix(ch_samples_without_reads.flatMap { meta -> [
            tuple(meta, 'preprocessing', [reads_left: 0], []),
            tuple(meta, 'genome_kmer_based_calling', [status: 'skipped_no_reads'], []),
            tuple(meta, 'genome_mapping_based_calling', [status: 'skipped_no_reads'], []),
        ] })
        .mix(params.skip_genome_kmer_based_calling
            ? ch_callable.map { meta, reads -> tuple(meta, 'genome_kmer_based_calling', [status: 'skipped'], []) }
            : ch_kmer_calls.map { meta, table -> tuple(meta, 'genome_kmer_based_calling', [calls_table: true], table) })
        .mix(params.skip_genome_mapping_based_calling
            ? ch_callable.map { meta, reads -> tuple(meta, 'genome_mapping_based_calling', [status: 'skipped'], []) }
            : ch_mapping_status.map { meta, status -> tuple(meta, 'genome_mapping_based_calling', [status: status], []) })
    WRITE_SAMPLE_STATUS(ch_status_rows, run_id)

    // ---- One MultiQC report per run for all read QC ----
    MULTIQC(
        ch_step_logs.map { sample_id, log_file -> log_file }.collect(),
        BUILD_PREPROCESSING_READ_COUNTS.out.table.map { meta, table -> table }.collect(),
        file("${projectDir}/assets/multiqc_config.yaml"),
        run_id,
    )

    emit:
    preprocessing_read_counts = BUILD_PREPROCESSING_READ_COUNTS.out.table
    preprocessing_logs = ch_step_logs.map { sample_id, log_file -> log_file }
    genome_kmer_based_calls = ch_kmer_calls
    genome_mapping_based_calls = ch_mapping_calls
    esviritu_reports = params.skip_genome_mapping_based_calling ? channel.empty() : ESVIRITU.out.reports
    esviritu_intermediates = params.skip_genome_mapping_based_calling ? channel.empty() : ESVIRITU.out.intermediates
    sample_status = WRITE_SAMPLE_STATUS.out.table
    qc_report = MULTIQC.out.report.mix(MULTIQC.out.data)
}
