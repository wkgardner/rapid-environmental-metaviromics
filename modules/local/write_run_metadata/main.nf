process WRITE_RUN_METADATA {
    // <run_id>.run_metadata.yaml: run identity, per-step tools and versions, every effective param.
    label 'process_single'

    input:
    path versions, stageAs: 'versions/*.yml'
    val run_id
    val params_json
    val start_time

    output:
    path "${run_id}.run_metadata.yaml", emit: metadata

    script:
    """
    cat > effective_params.json <<'END_PARAMS'
    ${params_json}
    END_PARAMS
    write_run_metadata.py \\
        --run-id ${run_id} --pipeline-version '${workflow.manifest.version}' --start-time '${start_time}' \\
        --params-json effective_params.json --versions versions/*.yml --out-dir .
    """
}
