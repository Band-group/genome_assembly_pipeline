process CREATE_BAM {
    tag "$meta.sample_id"
    //label "process_single"
    //publishDir "${projectDir}/results/03_BAM_files/", mode: "copy"

    input:
        tuple val(meta), path(assembly_fasta)

    output:
        tuple val(meta), path("${meta.sample_id}_sorted.bam"), emit: bam
        tuple val(meta), path("*.bai"), emit: bai

    script:
        """
        minimap2 -x asm5 -a ${meta.reference} ${assembly_fasta} > ${meta.sample_id}.sam
		samtools sort -o ${meta.sample_id}_sorted.bam ${meta.sample_id}.sam
		samtools index ${meta.sample_id}_sorted.bam
        """
}