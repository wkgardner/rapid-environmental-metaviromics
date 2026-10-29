process REFERENCE_FILTER_SRA {
    // NVD v3.5's SRA path: stream the run from `sracha get` straight into the same Deacon
    // reference-filter pass (deplete or extract, as REFERENCE_FILTER), so the decoded FASTQ
    // never lands on disk. Paired runs come out
    // interleaved, like local paired input. The run's layout is read by column name from
    // `sracha info` (sracha 0.8 added columns that NVD's fixed 8-column check rejects).
    tag "${meta.id}, ${meta.srr}"
    label 'process_medium'
    maxForks params.max_concurrent_downloads

    input:
    tuple val(meta), path(filter_index)
    val filter_mode

    output:
    tuple val(meta), path("${meta.id}.reference_filtered.fastq.gz"), env('READS_POST_REFERENCE_FILTER'), env('READ_MODE'), emit: reads
    tuple val(meta), path("${meta.id}.reference_filter.deacon.json"), emit: summary
    path 'versions.yml', topic: versions

    script:
    if (!(filter_mode in ['deplete', 'extract'])) {
        error "reference_filter_mode must be 'deplete' or 'extract', got '${filter_mode}'"
    }
    def deplete_argument = filter_mode == 'deplete' ? '--deplete' : ''
    def summary = "${meta.id}.reference_filter.deacon.json"
    def complexity_arg = params.reference_filter_complexity_threshold != null
        ? "--complexity-threshold ${params.reference_filter_complexity_threshold}" : ''
    def sracha_threads = Math.max(1, (task.cpus as int).intdiv(2))
    def deacon_threads = Math.max(1, (task.cpus as int) - sracha_threads)
    """
    set -euo pipefail

    sracha info --format tsv '${meta.srr}' > sracha_info.tsv
    layout=\$(awk -F'\\t' 'NR == 1 { for (i = 1; i <= NF; i++) column[\$i] = i } NR == 2 { print \$column["layout"] "/" \$column["nreads"] }' sracha_info.tsv)

    deacon_args=()
    case "\${layout}" in
        SINGLE/1) READ_MODE=single ;;
        PAIRED/2) READ_MODE=paired; deacon_args+=(--interleaved ${params.check_pairs ? '--check-pairs' : ''}) ;;
        *) echo "Unsupported SRA layout for ${meta.srr}: \${layout}" >&2; exit 1 ;;
    esac

    sracha get --output-dir . --stdout --split interleaved --threads ${sracha_threads} \\
        --no-progress --yes '${meta.srr}' \\
    | deacon filter ${deplete_argument} "\${deacon_args[@]}" \\
        --abs-threshold ${params.reference_filter_abs_threshold} \\
        --rel-threshold ${params.reference_filter_rel_threshold} ${complexity_arg} \\
        --threads ${deacon_threads} --summary ${summary} \\
        --output ${meta.id}.reference_filtered.fastq.gz \\
        ${filter_index} -

    READS_POST_REFERENCE_FILTER=\$(jq -er '.seqs_out' ${summary})
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deacon: "\$(deacon --version | sed 's/^deacon //')"
    END_VERSIONS
    """
}
