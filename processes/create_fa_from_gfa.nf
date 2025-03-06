process CREATE_FA_FROM_GFA {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/02_FASTA/", mode: "copy"

    container "oras://community.wave.seqera.io/library/gfatools_samtools:3251264bda7675f9" // singularity

    input: 
        tuple val(meta), path(assembly_graph)

    output:
        tuple val(meta), path("${meta.sample_id}.primary_contigs.fa.gz"), emit: contig_fa
    
    script:
        """
        gfatools gfa2fa ${assembly_graph} > ${meta.sample_id}.primary_contigs.fa
        bgzip ${meta.sample_id}.primary_contigs.fa
        """
}