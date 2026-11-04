#!/usr/bin/env nextflow
// Rapid wastewater metagenomics v1. Spec: docs/v1-spec.md. Usage:
//   nextflow run . -params-file my_run.yaml     (input, reference_filter_index, sylph_database,
//                                                vpd_lookup, esviritu_database, outdir, ...)

include { RAPID_WW           } from './workflows/rapid_ww'
include { WRITE_RUN_METADATA } from './modules/local/write_run_metadata'

workflow {
    main:
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
    genome_kmer_based_calls = RAPID_WW.out.genome_kmer_based_calls
    genome_mapping_based_calls = RAPID_WW.out.genome_mapping_based_calls
    esviritu_reports = RAPID_WW.out.esviritu_reports
    esviritu_intermediates = RAPID_WW.out.esviritu_intermediates
    sample_status = RAPID_WW.out.sample_status
    qc_report = RAPID_WW.out.qc_report
    run_metadata = WRITE_RUN_METADATA.out.metadata
}

// Published layout (spec "Outputs"): every table is one TSV per sample in <table>/.
// Files are copied, not symlinked (workflow.output.mode in nextflow.config).
output {
    preprocessing_read_counts { path '01_preprocessing/preprocessing_read_counts' }
    preprocessing_logs { path '01_preprocessing/logs' }
    genome_kmer_based_calls { path '02_genome_kmer_based_calling/genome_kmer_based_calls' }
    genome_mapping_based_calls { path '03_genome_mapping_based_calling/genome_mapping_based_calls' }
    // EsViritu's raw files land in <sample>/esviritu_output/ (the folder EsViritu writes them to).
    esviritu_reports { path { meta, report_files -> "03_genome_mapping_based_calling/esviritu_reports/${meta.id}" } }
    // With esviritu_keep on, the temp folder sits beside the raw files: <sample>/esviritu_output/<sample>_temp/.
    esviritu_intermediates { path { meta, temp_folder -> "03_genome_mapping_based_calling/esviritu_reports/${meta.id}" } }
    qc_report { path '04_qc_report' }
    sample_status { path '05_run_summary/sample_status' }
    run_metadata { path '05_run_summary' }
}
