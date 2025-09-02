process ALIGN_FA_TO_REF {
    tag "$meta.sample_id"
    label "process_single"

    container "oras://community.wave.seqera.io/library/minimap2_samtools:564e363a70ace2ad" // singularity

    input:
        tuple val(meta), path(assembly_fasta)

    output:
        tuple val(meta), path("${meta.sample_id}_sorted.bam"), emit: aligned_fa_bam
        tuple val(meta), path("*.bai"), emit: aligned_fa_bai
        path("versions.yml"), emit: versions

    script:
        def asm_preset = "${params.minimap_asm}"
        """ 
        minimap2 -x ${asm_preset} -a ${meta.reference} ${assembly_fasta} > ${meta.sample_id}.sam
		samtools sort -o ${meta.sample_id}_sorted.bam ${meta.sample_id}.sam
		samtools index ${meta.sample_id}_sorted.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            minimap2: \$(minimap2 --version 2>&1)
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}