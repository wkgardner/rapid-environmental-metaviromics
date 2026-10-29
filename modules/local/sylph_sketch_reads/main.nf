process SYLPH_SKETCH_READS {
    // Lifted from silly-fast-bfx 151c676 (MIT) SKETCH_READS_WITH_SYLPH: reads are sketched as one
    // single-end stream, which suits interleaved pairs and ONT reads alike (the spec's sylph -r).
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.sylsp"), emit: sketch
    path 'versions.yml', topic: versions

    script:
    """
    set -o pipefail
    zcat -f ${reads} \\
    | sylph sketch \\
        --reads /dev/stdin \\
        --sample-names ${meta.id} \\
        -k ${params.sylph_sample_kmer_length} \\
        -c ${params.sylph_sample_compression} \\
        -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: "\$(sylph --version | sed 's/^sylph //')"
    END_VERSIONS
    """
}
