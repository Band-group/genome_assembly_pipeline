process GET_CHROMOSOME_LENGTHS {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/01_FASTQ_QC_REPORTS/coverage/chromosome_lengths"

    input:
        tuple val(meta), path(bams)

    output:
        tuple val(meta), path(bams), path("chromosome_lengths.txt"), emit: chr_lengths

    script:
        """
        awk 'BEGIN {OFS = "\t"} {sum[\$1]+=\$3-\$2} END {for (chrom in sum) print chrom, sum[chrom]}' ${meta.regions_bed} > chromosome_lengths.txt
        """
}