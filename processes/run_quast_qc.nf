// Process to run Quast QC on oriented contigs

process RUN_QUAST_QC {
    tag "$meta.sample_id"
    label "process_low"
    publishDir "${params.outdir}/results/06_QUAST_QC/", mode: "copy"

    input:
        tuple val(meta), path(assembly_fasta)

    output:
        path("${meta.sample_id}_quast/*"), emit: quast_results

    script:
        """
        quast \\
            ${assembly_fasta} \\
            -o ${meta.sample_id}_quast \\
            -r ${meta.reference} \\
            -t ${task.cpus}
        """
}