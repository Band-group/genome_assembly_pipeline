// Process to run Quast QC on oriented contigs

process RUN_QUAST_QC {
    tag "$meta.sample_id"
    label "process_low"
    publishDir "${params.outdir}/results/06_QUAST_QC/", mode: "copy"

    input:
        tuple val(meta), path(assembly_fasta)

    output:
        path("quast_results/*"), emit: quast_results

    script:
        """
        quast \\
            ${assembly_fasta} \\
            -r ${meta.reference} \\
            -t ${task.cpus}
        """






}