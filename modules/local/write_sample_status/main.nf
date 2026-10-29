process WRITE_SAMPLE_STATUS {
    // One sample x step row of sample_status. status_source is a map with one of:
    // [status: 'failed'], [reads_left: 3584] or [calls_table: <path>] (passed as calls_table).
    tag "${meta.id}, ${step}"
    label 'process_single'

    input:
    tuple val(meta), val(step), val(status_source), path(calls_table, stageAs: 'calls/*')
    val run_id

    output:
    tuple val(meta), path("${run_id}_${meta.id}_${step}.sample_status.tsv"), emit: table

    script:
    def status_argument = status_source.status != null ? "--status ${status_source.status}"
        : status_source.reads_left != null ? "--reads-left ${status_source.reads_left}"
        : "--calls-table ${calls_table}"
    """
    write_sample_status.py --run-id ${run_id} --sample-id ${meta.id} --step ${step} ${status_argument} --out-dir .
    """
}
