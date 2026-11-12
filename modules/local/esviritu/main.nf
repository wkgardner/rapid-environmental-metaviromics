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
    // EsViritu's raw tables, consensus, log and params; not its HTML reports or their JavaScript folder.
    // Flattened to the task folder by the repackaging step below so they publish without a nested
    // esviritu_output/ level.
    tuple val(meta), path("*.{tsv,yaml,log,fasta,json}"), emit: reports, optional: true
    // The round-3 alignment: the one every downstream read-level question goes back to.
    tuple val(meta), path("${meta.id}.third.filt.sorted.bam*"), emit: final_bam, optional: true
    // Everything else EsViritu kept (rounds 1-2 and the clustering tables), when esviritu_keep is on.
    // Named for the sample alone so it publishes as esviritu_intermediates/<sample>/.
    tuple val(meta), path("${meta.id}", type: 'dir'), emit: intermediates, optional: true
    path 'versions.yml', topic: versions

    script:
    def read_format = meta.read_mode == 'paired' ? 'interleaved' : 'unpaired'
    def preset = meta.platform == 'illumina' ? params.esviritu_minimap2_preset_illumina : params.esviritu_minimap2_preset_nanopore
    // Kept files must stay in the task folder to be published; otherwise they can use fast local scratch.
    def temp_folder = params.esviritu_keep ? "esviritu_output/${meta.id}_temp" : '$TMPDIR/esviritu_temp'
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
        --keep ${params.esviritu_keep ? 'True' : 'False'} \\
        --temp "${temp_folder}" \\
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

    # Repackage EsViritu's layout for publishing. EsViritu writes its tables into esviritu_output/
    # and its working files into esviritu_output/<sample>_temp/. Published outputs want three flat,
    # self-describing groups instead, so move everything up to the task folder:
    #   <sample>.*                      tables, consensus, log, params        -> esviritu_results/
    #   <sample>.third.filt.sorted.bam  the round-3 alignment                 -> esviritu_final_bams/
    #   <sample>/                       rounds 1-2 and the clustering tables  -> esviritu_intermediates/
    if [ -d esviritu_output ]; then
        find esviritu_output -maxdepth 1 -type f -exec mv -f {} . \\;

        esviritu_temp_folder=esviritu_output/${meta.id}_temp
        if [ -d "\$esviritu_temp_folder" ]; then
            # The .bai may be absent if EsViritu never indexed; tolerate either.
            for final_alignment_file in "\$esviritu_temp_folder"/${meta.id}.third.filt.sorted.bam \\
                                        "\$esviritu_temp_folder"/${meta.id}.third.filt.sorted.bam.bai; do
                if [ -f "\$final_alignment_file" ]; then mv -f "\$final_alignment_file" .; fi
            done
            mv "\$esviritu_temp_folder" ${meta.id}
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esviritu: "\$(EsViritu --version 2>&1 | tail -1)"
    END_VERSIONS
    """
}
