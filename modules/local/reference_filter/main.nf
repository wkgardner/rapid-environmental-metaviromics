process REFERENCE_FILTER {
    // One Deacon pass against the reference-filter index. filter_mode "deplete" (default) drops
    // reads that match the index (host + rRNA); "extract" keeps only matching reads (e.g. VPD).
    // Paired input is written interleaved, so every later step sees one FASTQ per sample.
    // Lane-split or multi-file inputs stream through NVD's stream_fastqs_to_deacon.py.
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(read_files, stageAs: "reads??????/*")
    path filter_index
    val filter_mode

    output:
    tuple val(meta), path("${meta.id}.reference_filtered.fastq.gz"), env('READS_POST_REFERENCE_FILTER'), emit: reads
    tuple val(meta), path("${meta.id}.reference_filter.deacon.json"), emit: summary
    path 'versions.yml', topic: versions

    script:
    def files = read_files instanceof List ? read_files : [read_files]
    def r1_files = files.take(meta.r1_count as int)
    def r2_files = files.drop(meta.r1_count as int)
    if (!(filter_mode in ['deplete', 'extract'])) {
        error "reference_filter_mode must be 'deplete' or 'extract', got '${filter_mode}'"
    }
    def deplete_argument = filter_mode == 'deplete' ? '--deplete' : ''
    def summary = "${meta.id}.reference_filter.deacon.json"
    def output = "${meta.id}.reference_filtered.fastq.gz"
    def check_pairs_arg = params.check_pairs && meta.read_mode == 'paired' ? '--check-pairs' : ''
    def thresholds = "--abs-threshold ${params.reference_filter_abs_threshold} --rel-threshold ${params.reference_filter_rel_threshold}"
    def complexity_arg = params.reference_filter_complexity_threshold != null
        ? "--complexity-threshold ${params.reference_filter_complexity_threshold}" : ''
    def write_versions = """
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            deacon: "\$(deacon --version | sed 's/^deacon //')"
        END_VERSIONS""".stripIndent()
    def count_reads = "READS_POST_REFERENCE_FILTER=\$(jq -er '.seqs_out' ${summary})"
    if (files.size() == 1 || (meta.read_mode == 'paired' && r1_files.size() == 1))
        """
        deacon filter ${deplete_argument} ${check_pairs_arg} ${thresholds} ${complexity_arg} \\
            --threads ${task.cpus} --summary ${summary} --output ${output} \\
            ${filter_index} ${files.join(' ')}
        ${count_reads}
        ${write_versions}
        """
    else if (meta.read_mode == 'single')
        """
        printf '%s\\n' ${files.join(' ')} > reads.list
        stream_fastqs_to_deacon.py --sample-id ${meta.id} --index ${filter_index} ${deplete_argument} \\
            --reads-list reads.list ${thresholds} \\
            --threads ${task.cpus} --summary ${summary} --output ${output}
        ${count_reads}
        ${write_versions}
        """
    else
        """
        printf '%s\\n' ${r1_files.join(' ')} > r1.list
        printf '%s\\n' ${r2_files.join(' ')} > r2.list
        stream_fastqs_to_deacon.py --sample-id ${meta.id} --index ${filter_index} ${deplete_argument} ${check_pairs_arg} \\
            --r1-list r1.list --r2-list r2.list ${thresholds} \\
            --threads ${task.cpus} --summary ${summary} --output ${output}
        ${count_reads}
        ${write_versions}
        """
}
