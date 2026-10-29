process FETCH_SRA_READS {
    // Runs only when the reference filter is skipped: download an SRA run with sracha as one
    // FASTQ (paired runs interleaved), matching what REFERENCE_FILTER_SRA streams into Deacon.
    tag "${meta.id}, ${meta.srr}"
    label 'process_medium'
    maxForks params.max_concurrent_downloads

    input:
    val meta

    output:
    tuple val(meta), path("${meta.id}.prepared.fastq.gz"), env('READ_MODE'), emit: reads
    path 'versions.yml', topic: versions

    script:
    """
    set -euo pipefail

    sracha info --format tsv '${meta.srr}' > sracha_info.tsv
    layout=\$(awk -F'\\t' 'NR == 1 { for (i = 1; i <= NF; i++) column[\$i] = i } NR == 2 { print \$column["layout"] "/" \$column["nreads"] }' sracha_info.tsv)
    case "\${layout}" in
        SINGLE/1) READ_MODE=single ;;
        PAIRED/2) READ_MODE=paired ;;
        *) echo "Unsupported SRA layout for ${meta.srr}: \${layout}" >&2; exit 1 ;;
    esac

    sracha get --output-dir . --stdout --split interleaved --threads ${task.cpus} \\
        --no-progress --yes '${meta.srr}' \\
    | gzip -1 > ${meta.id}.prepared.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sracha: "\$(sracha --version | sed 's/^sracha //')"
    END_VERSIONS
    """
}
