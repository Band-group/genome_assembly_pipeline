process COMPUTE_COVERAGE_PER_BASE {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/03_CONTIG_CORRECTION/bedtools_genomecov/${meta.sample_id}", mode: "copy"
    
    container "oras://community.wave.seqera.io/library/bedtools:2.31.1--a120a7e98287539a"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("${meta.sample_id}_genomecov_coverage.txt"), emit: coverage_data
    
    script:
        """
        bedtools genomecov -ibam ${bam} -d > ${meta.sample_id}_genomecov_coverage.txt
        """
}