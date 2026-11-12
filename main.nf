#!/usr/bin/env nextflow
// Rapid wastewater metagenomics v1. Spec: docs/v1-spec.md. Usage:
//   nextflow run . -params-file my_run.yaml     (input, reference_filter_index, sylph_database,
//                                                esviritu_database, outdir, ...)

include { RAPID_WW           } from './workflows/rapid_ww'
include { WRITE_RUN_METADATA } from './modules/local/write_run_metadata'

workflow {
    main:
    // esviritu_keep decides how much of EsViritu's alignment work survives the run. Checked here
    // rather than deep in the process, so a typo fails in the first second instead of after the
    // last sample has mapped.
    def esviritu_keep_values = ['none', 'final', 'all']
    if (!(params.esviritu_keep in esviritu_keep_values)) {
        error "esviritu_keep must be one of ${esviritu_keep_values.join(', ')}, got '${params.esviritu_keep}'. " +
              "none: publish no alignments. final: the round-3 BAM, its index and its reference. " +
              "all: the above plus rounds 1-2 and the clustering tables."
    }
    // Every row of every table carries run_id; default rww_<YYYYMMDD>_<Nextflow runName>.
    def run_id = params.run_id ?: "rww_${workflow.start.format('yyyyMMdd')}_${workflow.runName}"
    RAPID_WW(file(params.input, checkIfExists: true), run_id)
    WRITE_RUN_METADATA(
        channel.topic('versions').collect(),
        run_id,
        groovy.json.JsonOutput.toJson(params),
        workflow.start.toString(),
    )

    publish:
    preprocessing_read_counts = RAPID_WW.out.preprocessing_read_counts
    preprocessing_logs = RAPID_WW.out.preprocessing_logs
    sylph_profiles = RAPID_WW.out.sylph_profiles
    esviritu_reports = RAPID_WW.out.esviritu_reports
    esviritu_final_bams = RAPID_WW.out.esviritu_final_bams
    esviritu_intermediates = RAPID_WW.out.esviritu_intermediates
    sample_status = RAPID_WW.out.sample_status
    qc_report = RAPID_WW.out.qc_report
    run_metadata = WRITE_RUN_METADATA.out.metadata
}

// Published layout (spec "Outputs"). Each caller publishes its own raw output; REMV adds
// only run-level bookkeeping (read counts, sample status, run metadata).
// Files are copied, not symlinked (workflow.output.mode in nextflow.config).
output {
    preprocessing_read_counts { path '01_preprocessing/preprocessing_read_counts' }
    preprocessing_logs { path '01_preprocessing/logs' }
    // sylph's own profile, one TSV per sample. REMV publishes what the caller wrote.
    sylph_profiles { path { meta, profile -> '02_genome_kmer_based_calling/sylph_results' } }
    // EsViritu's own tables, consensus, log and params, one folder per sample.
    esviritu_reports { path { meta, report_files -> "03_genome_mapping_based_calling/esviritu_results/${meta.id}" } }
    // The round-3 alignments, flat: filenames are already sample-prefixed, so one glob gets them all.
    esviritu_final_bams { path { meta, bam_files -> '03_genome_mapping_based_calling/esviritu_final_bams' } }
    // With esviritu_keep on, rounds 1-2 and the clustering tables, one folder per sample.
    esviritu_intermediates { path { meta, intermediates_folder -> '03_genome_mapping_based_calling/esviritu_intermediates' } }
    qc_report { path '04_qc_report' }
    sample_status { path '05_run_summary/sample_status' }
    run_metadata { path '05_run_summary' }
}
