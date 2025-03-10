process COMPUTE_COVERAGE_PER_BASE {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/02_CONTIG_ADJUSTMENT/bedtools_genomecov/${meta.sample_id}", mode: "copy"
    
    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("${meta.sample_id}_genomecov_coverage.txt"), emit: coverage_data
    
    script:
        """
        bedtools genomecov -ibam ${bam} -d > ${meta.sample_id}_genomecov_coverage.txt
        """
}