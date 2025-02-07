process COMPUTE_COVERAGE {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/01_FASTQ_QC_REPORTS/coverage/bedtools_genomecov_data", mode: "copy"
    
    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("${meta.sample_id}_genomecov_coverage.txt"), emit: coverage_data
    
    script:
        """
        bedtools genomecov -ibam ${bam} > ${meta.sample_id}_genomecov_coverage.txt
        """
}