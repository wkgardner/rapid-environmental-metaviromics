process BUILD_GENOME_KMER_BASED_CALLS {
    // One sample's genome_kmer_based_calls table from its sylph profile.
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(sylph_profile)
    path vpd_lookup
    val run_id

    output:
    tuple val(meta), path("${run_id}_${meta.id}.genome_kmer_based_calls.tsv"), emit: table

    script:
    """
    build_genome_kmer_based_calls.py \\
        --run-id ${run_id} --sample-id ${meta.id} \\
        --sylph-profile ${sylph_profile} --vpd-lookup ${vpd_lookup} --out-dir .
    """
}
