process MULTIQC {
    // One MultiQC report per run for all read QC: Deacon summaries and BBDuk logs through
    // MultiQC's own modules. Every sample's preprocessing_read_counts row fills the General
    // Statistics table (assets/multiqc_config.yaml), one row per sample with a column per step.
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
