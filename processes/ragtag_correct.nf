process RAGTAG_CORRECT {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/02_ASSEMBLIES/scaffolding/ragtag_correct/${meta.sample_id}"

    input:
        tuple val(meta), path(assembly_fasta) // output from orient_contigs

    output:
        tuple val(meta), path("ragtag.correct.fasta"), emit: rt_corrected_fa
        tuple val(meta), path("ragtag.correct*"), emit: rt_corrected_contigs

    script: // TODO: add ragtag parameters to config file
        """
        ragtag.py correct ${meta.reference} ${assembly_fasta} -o .
        """
}