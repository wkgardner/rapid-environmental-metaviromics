process REFERENCE_FILTER_CRAM {
    // Archived uCRAM path: stream `samtools fastq` straight into the same Deacon reference-filter
    // pass (deplete or extract, as REFERENCE_FILTER), so the expanded FASTQ never lands on disk.
    // Paired input comes out interleaved, like local paired input.
    //
    // A CRAM holds every lane and both mates of a sample in one file, and nothing in the
    // samplesheet says which, so the layout is read from the records' flags: any READ2 (0x80)
    // record means paired. required_fields=2 decodes flags only, skipping CRAM's sequence and
    // quality blocks, so that pass costs a fraction of the decode.
    //
    // samtools appends /1 and /2 to mate names, which is one of the two spellings
    // Deacon's --check-pairs accepts; -n would leave both mates with the same name and fail it.
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(cram), path(filter_index)
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
    def samtools_threads = Math.max(1, (task.cpus as int).intdiv(2))
    def deacon_threads = Math.max(1, (task.cpus as int) - samtools_threads)
    """
    set -euo pipefail

    read2_records=\$(samtools view -c -f 128 --input-fmt-option required_fields=2 ${cram})

    deacon_args=()
    if [ "\${read2_records}" -gt 0 ]; then
        READ_MODE=paired
        deacon_args+=(--interleaved ${params.check_pairs ? '--check-pairs' : ''})
    else
        READ_MODE=single
    fi

    samtools fastq -@ ${samtools_threads} ${cram} \\
    | deacon filter ${deplete_argument} "\${deacon_args[@]}" \\
        --abs-threshold ${params.reference_filter_abs_threshold} \\
        --rel-threshold ${params.reference_filter_rel_threshold} ${complexity_arg} \\
        --threads ${deacon_threads} --summary ${summary} \\
        --output ${meta.id}.reference_filtered.fastq.gz \\
        ${filter_index} -

    READS_POST_REFERENCE_FILTER=\$(jq -er '.seqs_out' ${summary})
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "\$(samtools --version | sed -n '1s/^samtools //p')"
        deacon: "\$(deacon --version | sed 's/^deacon //')"
    END_VERSIONS
    """
}
