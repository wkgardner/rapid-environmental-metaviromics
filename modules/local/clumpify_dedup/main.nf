process CLUMPIFY_DEDUP {
    // NVD's DEDUP_WITH_CLUMPIFY: remove duplicate reads (both mates must match for pairs).
    // subs=0 and allowns=f make it exact-match only, the same rule NVD uses.
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.deduplicated.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.dedup.clumpify.log"), emit: log
    path 'versions.yml', topic: versions

    script:
    def interleaved_argument = meta.read_mode == 'paired' ? 'int=t' : ''
    def java_heap_gigabytes = Math.max(1, (task.memory.toGiga() * 0.85) as int)
    """
    clumpify.sh -Xmx${java_heap_gigabytes}g \\
        in=${reads} \\
        out=${meta.id}.deduplicated.fastq.gz \\
        dedupe=t \\
        optical=f \\
        reorder=f \\
        subs=${params.dedup_substitutions} \\
        allowns=f \\
        ${interleaved_argument} \\
        usetmpdir=t tmpdir=\$TMPDIR \\
        threads=${task.cpus} \\
        -eoom \\
        2> >(tee ${meta.id}.dedup.clumpify.log >&2)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: "\$(bbversion.sh 2>/dev/null | tail -1)"
    END_VERSIONS
    """
}
