// Concatenate read files into one gzipped FASTQ: gzip members join cleanly, plain files are compressed.
def joinReads(reads, output_name) {
    def all_gzipped = reads.every { read_file -> read_file.name.endsWith('.gz') }
    return all_gzipped ? "cat ${reads.join(' ')} > ${output_name}" : "cat ${reads.join(' ')} | gzip -1 > ${output_name}"
}

process PREPARE_READS {
    // Runs only when the reference filter is skipped. Deacon normally turns each sample's
    // lane-split R1/R2 files into one interleaved FASTQ; this does the same without filtering:
    // lanes are concatenated (gzip members join cleanly) and pairs interleaved with reformat.sh.
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(read_files, stageAs: "reads??????/*")

    output:
    tuple val(meta), path("${meta.id}.prepared.fastq.gz"), emit: reads
    path 'versions.yml', topic: versions

    script:
    def files = read_files instanceof List ? read_files : [read_files]
    def r1_files = files.take(meta.r1_count as int)
    def r2_files = files.drop(meta.r1_count as int)
    def java_heap_gigabytes = Math.max(1, (task.memory.toGiga() * 0.85) as int)
    def prepare_command = meta.read_mode == 'paired'
        ? """
        ${joinReads(r1_files, 'joined_R1.fastq.gz')}
        ${joinReads(r2_files, 'joined_R2.fastq.gz')}
        reformat.sh -Xmx${java_heap_gigabytes}g in1=joined_R1.fastq.gz in2=joined_R2.fastq.gz \\
            out=${meta.id}.prepared.fastq.gz threads=${task.cpus} -eoom
        rm -f joined_R1.fastq.gz joined_R2.fastq.gz
        """
        : joinReads(files, "${meta.id}.prepared.fastq.gz")
    """
    ${prepare_command}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: "\$(bbversion.sh 2>/dev/null | tail -1)"
    END_VERSIONS
    """
}
