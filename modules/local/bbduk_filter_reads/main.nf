process BBDUK_FILTER_READS {
    // NVD's FILTER_READS: drop reads by average quality and length, and low-complexity reads by
    // k-mer entropy. For pairs, a pair is dropped when either mate fails (removeifeitherbad).
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.read_filtered.fastq.gz"), env('READS_POST_QC'), emit: reads
    tuple val(meta), path("${meta.id}.read_filtering.bbduk.log"), emit: log

    script:
    def min_quality = meta.platform == 'illumina' ? params.min_read_quality_illumina : params.min_read_quality_nanopore
    def entropy_arguments = params.filter_low_complexity_reads
        ? "entropy=${params.min_read_entropy} entropywindow=${params.min_read_entropy_window} entropyk=${params.min_read_entropy_kmer}"
        : ''
    def pair_arguments = meta.read_mode == 'paired' ? 'int=t removeifeitherbad=t' : 'int=f'
    def java_heap_gigabytes = Math.max(1, (task.memory.toGiga() * 0.85) as int)
    """
    bbduk.sh -Xmx${java_heap_gigabytes}g \\
        in=${reads} \\
        out=${meta.id}.read_filtered.fastq.gz \\
        ${pair_arguments} \\
        minavgquality=${min_quality} \\
        minlength=${params.min_read_length} \\
        ${entropy_arguments} \\
        ordered=t \\
        threads=${task.cpus} \\
        -eoom \\
        2> >(tee ${meta.id}.read_filtering.bbduk.log >&2)

    READS_POST_QC=\$(( \$(zcat ${meta.id}.read_filtered.fastq.gz | wc -l) / 4 ))
    """
}
