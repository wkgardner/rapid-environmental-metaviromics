// "--flag <log>" for a step that ran, "--flag skipped" for one a param switched off, else nothing.
def stepArgument(logs, flag, suffix, switched_off) {
    def log = logs.find { log_file -> log_file.name.endsWith(suffix) }
    return switched_off ? "${flag} skipped" : (log ? "${flag} ${log}" : '')
}

process BUILD_PREPROCESSING_READ_COUNTS {
    // One sample's preprocessing_read_counts row from whichever step logs exist. A step that a
    // param switched off reads "skipped"; a step that did not run for this sample (no reads
    // left, or adapter trimming on ONT) is left empty.
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(step_logs, stageAs: 'step_logs/*')
    val run_id

    output:
    tuple val(meta), path("${run_id}_${meta.id}.preprocessing_read_counts.tsv"), emit: table

    script:
    def logs = step_logs instanceof List ? step_logs : [step_logs]
    """
    build_preprocessing_read_counts.py \\
        --run-id ${run_id} --sample-id ${meta.id} --out-dir . \\
        ${stepArgument(logs, '--reference-filter-summary', '.reference_filter.deacon.json', params.skip_reference_filter)} \\
        ${stepArgument(logs, '--adapter-trimming-log', '.adapter_trimming.bbduk.log', !params.trim_adapters)} \\
        ${stepArgument(logs, '--dedup-log', '.dedup.clumpify.log', !params.dedup)} \\
        ${stepArgument(logs, '--read-filtering-log', '.read_filtering.bbduk.log', params.skip_read_filtering)}
    """
}
