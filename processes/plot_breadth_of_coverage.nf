process PLOT_BREADTH_OF_COVERAGE {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/01_FASTQ_QC_REPORTS/deeptools_coverage_plot"



}