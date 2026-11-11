process DECODE_CRAM_READS {
    // Turn one archived unaligned CRAM back into FASTQ, so retrospective runs can read cold
    // storage directly. A CRAM holds every lane and both mates of a sample in one file, and
    // nothing in the samplesheet says which, so the pairing is read from the records' flags:
    // any READ2 (0x80) record means the sample is paired. Reading flags alone skips CRAM's
    // sequence and quality blocks, so the check is cheap next to the decode itself.
    // Read names are kept exactly as archived (-n), rather than having /1 and /2 appended.
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(cram)

    output:
    tuple val(meta), path("${meta.id}.cram_reads*.fastq.gz"), env('READ_MODE'), emit: reads
    path 'versions.yml', topic: versions

    script:
    """
    set -euo pipefail

    read2_records=\$(samtools view -c -f 128 --input-fmt-option required_fields=2 ${cram})
    if [ "\${read2_records}" -gt 0 ]; then
        READ_MODE=paired
        samtools fastq -@ ${task.cpus} -n \\
            -1 ${meta.id}.cram_reads_R1.fastq.gz \\
            -2 ${meta.id}.cram_reads_R2.fastq.gz \\
            -0 /dev/null -s /dev/null ${cram}
    else
        READ_MODE=single
        # Unpaired records land in -0, and records flagged as one mate of a pair whose mate is
        # absent land in -s. Both are single reads here, and gzip members concatenate cleanly.
        samtools fastq -@ ${task.cpus} -n \\
            -0 unpaired.fastq.gz -s singleton.fastq.gz ${cram}
        cat unpaired.fastq.gz singleton.fastq.gz > ${meta.id}.cram_reads.fastq.gz
        rm -f unpaired.fastq.gz singleton.fastq.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "\$(samtools --version | sed -n '1s/^samtools //p')"
    END_VERSIONS
    """
}
