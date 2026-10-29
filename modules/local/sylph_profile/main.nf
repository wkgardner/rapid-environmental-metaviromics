process SYLPH_PROFILE {
    // Lifted from silly-fast-bfx 151c676 (MIT) PROFILE_READS_WITH_SYLPH, plus the hidden
    // -R redundancy threshold. Every output row is a genome k-mer-based call.
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(sketch)
    path database

    output:
    tuple val(meta), path("${meta.id}.sylph_profile.tsv"), emit: profile
    path 'versions.yml', topic: versions

    script:
    """
    sylph profile \\
        ${database} \\
        ${sketch} \\
        --minimum-ani ${params.sylph_minimum_ani} \\
        --min-count-correct ${params.sylph_min_count_correct} \\
        --min-number-kmers ${params.sylph_min_reference_kmers} \\
        --min-contain ${params.sylph_min_matching_kmers} \\
        -R ${params.sylph_redundancy_threshold} \\
        -t ${task.cpus} \\
        --output-file ${meta.id}.sylph_profile.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: "\$(sylph --version | sed 's/^sylph //')"
    END_VERSIONS
    """
}
