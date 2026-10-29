process RESOLVE_READ_INPUTS {
    // NVD's samplesheet rules (bin/nvd_read_inputs.py): one JSON record per sample.
    tag "${samplesheet.name}"
    label 'process_single'

    input:
    path samplesheet

    output:
    path "resolved_reads.jsonl", emit: jsonl

    script:
    """
    resolve_read_inputs.py \\
        --samplesheet ${samplesheet} \\
        --output-jsonl resolved_reads.jsonl
    """
}
