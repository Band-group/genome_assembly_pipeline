process COMPUTE_COVERAGE {
    tag "$meta.sample_id"
    label "process_low"
    publishDir "${params.outdir}/results/01_FASTQ_QC_REPORTS/bedtools_genomecov_data"
    
    input:
        tuple val(meta), path(bam), path(chr_lengths)

    output:
        tuple val(meta), path("${meta.sample_id}_coverage.txt"), emit: coverage_data
    
    script:
        if(chr_lengths){
            """
            bedtools bamtobed -i ${bam} \
            | bedtools genomecov -g ${chr_lengths} > ${meta.sample_id}_coverage.txt
            """
        } else {
            """
            bedtools genomecov -ibam ${bam} > ${meta.sample_id}_coverage.txt
            """
        }
}