process ESVIRITU {
    // Genome mapping-based calling. EsViritu's exit codes are unreliable, so the wrapper decides
    // the sample's status from what EsViritu left behind (spec):
    //   tables written                            -> ok
    //   exit 1 after "No reads aligned"           -> no_detections (a real negative)
    //   anything else (crash, exit 0, no tables)  -> failed
    // The task itself always succeeds so one sample cannot stop the run; the status goes to
    // sample_status instead.
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path database

    output:
    tuple val(meta), path("${meta.id}.detected_virus.assembly_summary.tsv"), path("${meta.id}.detected_virus.info.tsv"), env('MAPPING_STATUS'), emit: tables, optional: true
    tuple val(meta), env('MAPPING_STATUS'), emit: status
    tuple val(meta), path("${meta.id}_esviritu.log"), emit: log, optional: true
    path "esviritu_output/*", emit: reports, optional: true
    path 'versions.yml', topic: versions

    script:
    def read_format = meta.read_mode == 'paired' ? 'interleaved' : 'unpaired'
    def preset = meta.platform == 'illumina' ? params.esviritu_minimap2_preset_illumina : params.esviritu_minimap2_preset_nanopore
    """
    set +e
    EsViritu \\
        --reads ${reads} \\
        --read_format ${read_format} \\
        --sample ${meta.id} \\
        --output_dir esviritu_output \\
        --db ${database} \\
        --cpu ${task.cpus} \\
        --qual False \\
        --filter_seqs False \\
        --dedup ${params.esviritu_dedup ? 'True' : 'False'} \\
        --keep False \\
        --temp "\$TMPDIR/esviritu_temp" \\
        --minimap2-preset ${preset} \\
        --minimap2-K ${params.esviritu_minimap2_batch} \\
        --min-align-length ${params.esviritu_min_align_length} \\
        --species-threshold ${params.esviritu_species_threshold} \\
        --subspecies-threshold ${params.esviritu_subspecies_threshold}
    esviritu_exit_code=\$?
    set -e

    log_file=esviritu_output/${meta.id}_esviritu.log
    [ -f "\$log_file" ] && cp "\$log_file" ${meta.id}_esviritu.log
    if [ -f esviritu_output/${meta.id}.detected_virus.assembly_summary.tsv ] && [ -f esviritu_output/${meta.id}.detected_virus.info.tsv ]; then
        MAPPING_STATUS=ok
        cp esviritu_output/${meta.id}.detected_virus.assembly_summary.tsv esviritu_output/${meta.id}.detected_virus.info.tsv .
    elif [ "\$esviritu_exit_code" -eq 1 ] && grep -q 'No reads aligned' "\$log_file" 2>/dev/null; then
        MAPPING_STATUS=no_detections
    else
        MAPPING_STATUS=failed
        echo "EsViritu failed for ${meta.id} (exit \$esviritu_exit_code, no tables)" >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esviritu: "\$(EsViritu --version 2>&1 | tail -1)"
    END_VERSIONS
    """
}
