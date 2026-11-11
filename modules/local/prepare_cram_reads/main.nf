process PREPARE_CRAM_READS {
    // Runs only when the reference filter is skipped: decode an archived uCRAM to the one
    // FASTQ per sample every later step expects (paired reads interleaved), matching what
    // REFERENCE_FILTER_CRAM streams into Deacon. The layout is read from the records' flags,
    // as there too.
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(cram)

    output:
    tuple val(meta), path("${meta.id}.prepared.fastq.gz"), env('READ_MODE'), emit: reads
    path 'versions.yml', topic: versions

    script:
    def samtools_threads = Math.max(1, (task.cpus as int) - 1)
    """
    set -euo pipefail

    read2_records=\$(samtools view -c -f 128 --input-fmt-option required_fields=2 ${cram})
    if [ "\${read2_records}" -gt 0 ]; then
        READ_MODE=paired
    else
        READ_MODE=single
    fi

    samtools fastq -@ ${samtools_threads} ${cram} | gzip -1 > ${meta.id}.prepared.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "\$(samtools --version | sed -n '1s/^samtools //p')"
    END_VERSIONS
    """
}
