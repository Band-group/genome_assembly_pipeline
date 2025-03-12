process COMPUTE_COV_ZSCORES {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/03_CONTIG_CORRECTION/bedtools_genomecov_zscores/${meta.sample_id}", mode: "copy"

    container "oras://community.wave.seqera.io/library/r-argparse_r-dplyr_r-ggplot2:4c31e454ced817ad" // singularity

    input:
        tuple val(meta), path(coverage_data)

    output:
        path("high_zscore_regions.bed"), emit: zscore_bed
    
    script:
        """
        compute_cov_zscores.r --coverage_file ${coverage_data} --zscore_threshold 3.0 --outdir ./
        """
}