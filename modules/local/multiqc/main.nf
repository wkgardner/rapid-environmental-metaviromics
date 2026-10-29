process MULTIQC {
    // One MultiQC report per run for all read QC: Deacon summaries and BBDuk logs through
    // MultiQC's own modules, plus every sample's preprocessing_read_counts row as a custom
    // table (MultiQC has no Clumpify module, so this is where duplicate counts appear).
    label 'process_single'

    input:
    path qc_logs, stageAs: 'qc_logs/*'
    path read_count_tables, stageAs: 'read_counts/*'
    path multiqc_config
    val run_id

    output:
    path "${run_id}.multiqc_report.html", emit: report
    path "${run_id}.multiqc_report_data", emit: data
    path 'versions.yml', topic: versions

    script:
    """
    {
        printf '# id: preprocessing_read_counts\\n'
        printf '# section_name: Preprocessing read counts\\n'
        printf '# description: One row per sample, as in the preprocessing_read_counts table.\\n'
        printf '# plot_type: table\\n'
        head -n 1 "\$(ls read_counts/* | head -n 1)" | cut -f2-
        for table in read_counts/*; do tail -n +2 "\$table" | cut -f2-; done
    } > preprocessing_read_counts_mqc.tsv

    multiqc --no-ai --force --config ${multiqc_config} \\
        --filename ${run_id}.multiqc_report.html \\
        qc_logs preprocessing_read_counts_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: "\$(multiqc --version | sed 's/^multiqc, version //')"
    END_VERSIONS
    """
}
