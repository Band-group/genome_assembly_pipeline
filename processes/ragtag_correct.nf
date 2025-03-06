process RAGTAG_CORRECT {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/02_ASSEMBLIES/scaffolding/ragtag_correct/${meta.sample_id}", mode: "copy"

    container "oras://community.wave.seqera.io/library/minimap2_ragtag:5ad34249839dbbbc"

    input:
        tuple val(meta), path(assembly_fasta) // output from orient_contigs
        tuple val(meta2), path(fastq)         // raw fastq reads

    output:
        tuple val(meta), path("ragtag.correct.fasta"), emit: rt_corrected_fa
        tuple val(meta), path("ragtag.correct*"), emit: rt_corrected_contigs
        path("versions.yml"), emit: versions

    script: // TODO: add ragtag parameters to config file
        def VERSION = '2.1.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
        def mm2 = "${params.minimap_asm}"
        """
        gzip -cd ${fastq} | sed -n '1~4s/^@/>/p;2~4p' > ${meta.sample_id}.fasta
        ragtag.py correct ${meta.reference} ${assembly_fasta} -R ${meta.sample_id}.fasta --aligner minimap2 --mm2-params ${mm2} -T corr -o .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ragtag: ${VERSION}
        END_VERSIONS
        """
}