process PLOT_COVERAGE {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/01_FASTQ_QC_REPORTS/coverage/plot"

    input:
        val(meta)
        path(coverage_data)
        val(ref_length)

    output:
        path("coverage_plot.pdf"), emit: coverage_plot
    
    script:
        if ( params.use_region ) {
            """
            plot_coverage.r --outdir ./ --sample_names ${meta.sample_id} --coverage_files ${coverage_data} --regions_bed ${meta.regions_bed[0]} --ref_length ${ref_length}
            """
        } else {
            """
            plot_coverage.r --outdir ./ --sample_names ${meta.sample_id} --coverage_files ${coverage_data}
            """
        }
        
}