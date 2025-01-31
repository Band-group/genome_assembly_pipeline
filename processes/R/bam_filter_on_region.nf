process BAM_FILTER_ON_REGION {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/01_FASTQ_QC_REPORTS/coverage/reads2ref_region_restricted", mode: "copy"

    input:
        tuple val(meta), path(reads_bam) // reads_bam is the downsampled bam file

    output:
        tuple val(meta), path("*regionRestricted.*"), emit: reads_bam_filtered

    script:
        def regions = "--regions_bed ${params.regions_bed}"
        """
        Rscript bam_filter_on_region.r --outdir . $regions --bams ${reads_bam}
        """
}