process PLOT_BREADTH_OF_COVERAGE {
    label "process_medium"
    publishDir "${params.outdir}/results/01_FASTQ_QC_REPORTS/deeptools_coverage_plot"

    input:
        val(bamTuplesList) // collected list of tuples => [ [meta, bamPath1], [meta, bamPath2], ... ]
        path(regions_bed)
    
    output:
        path("regionRestricted_coverage.png"), emit: coverage_png

    script:
        def bamPaths = bamTuplesList.collect{ it[1] }.join(' ')
        """
        plotCoverage -b ${bamPaths} --BED ${regions_bed} -p ${task.cpus} -o regionRestricted_coverage.png
        """
}