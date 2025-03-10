process BAM_FILTER_ON_REGION {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/01_FASTQ_QC_REPORTS/coverage/reads2ref_region_restricted/${meta.sample_id}", mode: "copy"

    container "oras://community.wave.seqera.io/library/bioconductor-genomicalignments_bioconductor-genomicranges_bioconductor-rsamtools_bioconductor-rtracklayer_r-argparse:b3bf752c05aefd1d"

    input:
        tuple val(meta), path(bam) // reads_bam is the downsampled bam file

    output:
        tuple val(meta), path("*regionRestricted*.bam"), emit: aligned_reads_bam_filtered
        tuple val(meta), path("*regionRestricted*.bai"), emit: aligned_reads_bai_filtered

    script:
        """
        bam_filter_on_region.r --outdir . --regions_bed ${meta.regions_bed} --bams ${bam}
        """
}