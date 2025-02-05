process PLOT_COVERAGE {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/01_FASTQ_QC_REPORTS/coverage/plot"

    input:
        val(meta)
        path(coverage_data)

    output:
        path("coverage_plot.pdf"), emit: coverage_plot
    
    script:
        """
        coverage_plot.r --outdir ./ --sample_names ${meta.sample_id} --coverage_files ${coverage_data}
        """
}