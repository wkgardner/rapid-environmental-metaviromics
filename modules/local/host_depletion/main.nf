process HOST_DEPLETION {
    // One Deacon `filter --deplete` pass against the unioned host + rRNA index.
    // Paired input is written interleaved, so every later step sees one FASTQ per sample.
    // Lane-split or multi-file inputs stream through NVD's stream_fastqs_to_deacon.py.
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(read_files, stageAs: "reads??????/*")
    path depletion_index

    output:
    tuple val(meta), path("${meta.id}.host_depleted.fastq.gz"), env('READS_POST_DEPLETION'), emit: reads
    tuple val(meta), path("${meta.id}.host_depletion.deacon.json"), emit: summary

    script:
    def files = read_files instanceof List ? read_files : [read_files]
    def r1_files = files.take(meta.r1_count as int)
    def r2_files = files.drop(meta.r1_count as int)
    def summary = "${meta.id}.host_depletion.deacon.json"
    def output = "${meta.id}.host_depleted.fastq.gz"
    def check_pairs_arg = params.check_pairs && meta.read_mode == 'paired' ? '--check-pairs' : ''
    def thresholds = "--abs-threshold ${params.host_depletion_abs_threshold} --rel-threshold ${params.host_depletion_rel_threshold}"
    def complexity_arg = params.host_depletion_complexity_threshold != null
        ? "--complexity-threshold ${params.host_depletion_complexity_threshold}" : ''
    def count_reads = "READS_POST_DEPLETION=\$(jq -er '.seqs_out' ${summary})"
    if (files.size() == 1 || (meta.read_mode == 'paired' && r1_files.size() == 1))
        """
        deacon filter --deplete ${check_pairs_arg} ${thresholds} ${complexity_arg} \\
            --threads ${task.cpus} --summary ${summary} --output ${output} \\
            ${depletion_index} ${files.join(' ')}
        ${count_reads}
        """
    else if (meta.read_mode == 'single')
        """
        printf '%s\\n' ${files.join(' ')} > reads.list
        stream_fastqs_to_deacon.py --sample-id ${meta.id} --index ${depletion_index} --deplete \\
            --reads-list reads.list ${thresholds} \\
            --threads ${task.cpus} --summary ${summary} --output ${output}
        ${count_reads}
        """
    else
        """
        printf '%s\\n' ${r1_files.join(' ')} > r1.list
        printf '%s\\n' ${r2_files.join(' ')} > r2.list
        stream_fastqs_to_deacon.py --sample-id ${meta.id} --index ${depletion_index} --deplete ${check_pairs_arg} \\
            --r1-list r1.list --r2-list r2.list ${thresholds} \\
            --threads ${task.cpus} --summary ${summary} --output ${output}
        ${count_reads}
        """
}
