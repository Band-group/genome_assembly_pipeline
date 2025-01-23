process ALIGN_READS_TO_REF {
    tag "$meta.sample_id"
    label "process_low"
    publishDir "${params.outdir}/results/02_READS_ALIGNED_TO_REF/", mode: "copy"

    input:
        tuple val(meta), path(fastq)

    output:
        tuple val(meta), path("${meta.sample_id}_reads2ref_sorted.bam"), emit: aligned_reads_bam
        tuple val(meta), path("*.bai"), emit: aligned_reads_bai

    script:
        """
        minimap2 -ax map-hifi -t ${task.cpus} ${fastq} ${meta.reference} > ${meta.sample_id}_reads2ref.sam
        samtools sort -o ${meta.sample_id}_reads2ref_sorted.bam ${meta.sample_id}_reads2ref.sam
		samtools index ${meta.sample_id}_reads2ref_sorted.bam
        """
}