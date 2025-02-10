process RAGTAG_SCAFFOLD {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/02_ASSEMBLIES/scaffolding/ragtag_scaffold/${meta.sample_id}"

    input:
        tuple val(meta), path(corrected_fasta) // output from orient_contigs

    output:
        tuple val(meta), path("ragtag.scaffold.fasta"), emit: rt_scaffolds_fa
        tuple val(meta), path("ragtag.scaffold*"), emit: rt_scaffolds

    script: // TODO: add ragtag parameters to config file
        """
        ragtag.py scaffold ${meta.reference} ${corrected_fasta} -o .
        """

}