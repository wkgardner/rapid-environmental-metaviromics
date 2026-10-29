process BUILD_GENOME_MAPPING_BASED_CALLS {
    // One sample's genome_mapping_based_calls table from EsViritu's tables (header only when
    // EsViritu reported no_detections or failed).
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), val(mapping_status), path(esviritu_tables, stageAs: 'esviritu/*')
    path vpd_lookup
    val run_id

    output:
    tuple val(meta), path("${run_id}_${meta.id}.genome_mapping_based_calls.tsv"), emit: table

    script:
    def tables = esviritu_tables instanceof List ? esviritu_tables : [esviritu_tables]
    def assembly_summary = tables.find { it.name.endsWith('.detected_virus.assembly_summary.tsv') }
    def info = tables.find { it.name.endsWith('.detected_virus.info.tsv') }
    def table_arguments = mapping_status == 'ok' ? "--assembly-summary ${assembly_summary} --info ${info}" : ''
    """
    build_genome_mapping_based_calls.py \\
        --run-id ${run_id} --sample-id ${meta.id} --mapping-status ${mapping_status} ${table_arguments} \\
        --vpd-lookup ${vpd_lookup} \\
        --min-read-count ${params.mapping_based_call_min_read_count} \\
        --species-threshold ${params.esviritu_species_threshold} \\
        --subspecies-threshold ${params.esviritu_subspecies_threshold} \\
        --normalization-genus '${params.normalization_genus}' \\
        --out-dir .
    """
}
