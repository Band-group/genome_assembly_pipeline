process COMPUTE_COVERAGE {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/01_FASTQ_QC_REPORTS/coverage/bedtools_genomecov/${meta.sample_id}", mode: "copy"
    
    container "oras://community.wave.seqera.io/library/bedtools:2.31.1--a120a7e98287539a"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("${meta.sample_id}_genomecov_coverage.txt"), emit: coverage_data
        path("versions.yml"), emit: versions
    
    script:
        """
        bedtools genomecov -ibam ${bam} > ${meta.sample_id}_genomecov_coverage.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version 2>&1)
        END_VERSIONS
        """
}