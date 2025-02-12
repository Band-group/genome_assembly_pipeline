process RAGTAG_CORRECT {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/02_ASSEMBLIES/scaffolding/ragtag_correct/${meta.sample_id}", mode: "copy"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ragtag:2.1.0--pyhb7b1952_0'
        : 'biocontainers/ragtag:2.1.0--pyhb7b1952_0'}"

    input:
        tuple val(meta), path(assembly_fasta) // output from orient_contigs

    output:
        tuple val(meta), path("ragtag.correct.fasta"), emit: rt_corrected_fa
        tuple val(meta), path("ragtag.correct*"), emit: rt_corrected_contigs
        path("versions.yml"), emit: versions

    script: // TODO: add ragtag parameters to config file
        def VERSION = '2.1.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
        """
        ragtag.py correct ${meta.reference} ${assembly_fasta} -o .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ragtag: ${VERSION}
        END_VERSIONS
        """
}