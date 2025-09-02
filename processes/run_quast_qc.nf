// Process to run Quast QC on oriented contigs

process RUN_QUAST_QC {
    tag "$meta.sample_id"
    label "process_low"
    publishDir "${params.outdir}/results/04_FINAL_QC/quast_qc/${meta.sample_id}", mode: "copy"

    container "oras://community.wave.seqera.io/library/quast:5.3.0--bfd4c029fde7e696" // singularity

    input:
        tuple val(meta), path(assembly_fasta)

    output:
        path("${meta.sample_id}_quast/*"), emit: quast_results
        path("versions.yml"), emit: versions

    script:
        """
        quast \\
            ${assembly_fasta} \\
            -o ${meta.sample_id}_quast \\
            -r ${meta.reference} \\
            -t ${task.cpus}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            quast: \$(quast --version 2>&1)
        END_VERSIONS
        """
}