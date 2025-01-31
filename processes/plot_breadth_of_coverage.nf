process PLOT_BREADTH_OF_COVERAGE {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/01_FASTQ_QC_REPORTS/deeptools_coverage_plot"

    input:
        val(bamTuplesList) // collected list of tuples => [ [meta, bamPath1], [meta, bamPath2], ... ]
        path(regions_bed)
    
    output:
        path("regionRestricted_coverage.png"), emit: coverage_png

    script:
        def bamPaths = bamTuplesList*.second.join(' ') // Using Groovy's spread operator to collect the second element of each tuple
        """
        plotCoverage -b ${bamPaths} --BED ${regions_bed} -o regionRestricted_coverage.png
        """
}