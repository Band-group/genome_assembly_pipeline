process ALIGN_FA_TO_REF {
    tag "$meta.sample_id"
    label "process_single"

    container "oras://community.wave.seqera.io/library/minimap2_samtools:564e363a70ace2ad" // singularity

    input:
        tuple val(meta), path(assembly_fasta)

    output:
        tuple val(meta), path("${meta.sample_id}_sorted.bam"), emit: aligned_fa_bam
        tuple val(meta), path("*.bai"), emit: aligned_fa_bai

    script:
        def asm_preset = "${params.minimap_asm}"
        """ 
        minimap2 -x ${asm_preset} -a ${meta.reference} ${assembly_fasta} > ${meta.sample_id}.sam
		samtools sort -o ${meta.sample_id}_sorted.bam ${meta.sample_id}.sam
		samtools index ${meta.sample_id}_sorted.bam
        """
}