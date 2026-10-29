process BBDUK_TRIM_ADAPTERS {
    // NVD's TRIM_ADAPTERS: right-trim Illumina adapters by k-mer match against BBDuk's bundled
    // adapter set; for pairs, also trim by pair overlap (tbo, tpe). Illumina reads only.
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.adapter_trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.adapter_trimming.bbduk.log"), emit: log

    script:
    def pair_arguments = meta.read_mode == 'paired' ? 'int=t tpe tbo removeifeitherbad=t' : 'int=f'
    def java_heap_gigabytes = Math.max(1, (task.memory.toGiga() * 0.85) as int)
    """
    bbduk.sh -Xmx${java_heap_gigabytes}g \\
        in=${reads} \\
        out=${meta.id}.adapter_trimmed.fastq.gz \\
        ref=adapters \\
        ktrim=r \\
        k=${params.adapter_trim_kmer_length} \\
        mink=${params.adapter_trim_min_kmer_length} \\
        hdist=${params.adapter_trim_hamming_distance} \\
        ${pair_arguments} \\
        threads=${task.cpus} \\
        -eoom \\
        2> >(tee ${meta.id}.adapter_trimming.bbduk.log >&2)
    """
}
