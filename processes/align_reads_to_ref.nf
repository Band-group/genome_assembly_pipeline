process ALIGN_READS_TO_REF {
    tag "$meta.sample_id"
    label "process_low"
    publishDir "${params.outdir}/results/01_FASTQ_QC_REPORTS/coverage/reads2ref", mode: "copy"

    container "oras://community.wave.seqera.io/library/minimap2_samtools:564e363a70ace2ad" // singularity

    input:
        tuple val(meta), path(fastq)

    output:
        tuple val(meta), path("${meta.sample_id}_reads2ref_sorted.bam"), emit: aligned_reads_bam
        tuple val(meta), path("*.bai"), emit: aligned_reads_bai

    script:
        """
        minimap2 -ax map-hifi -t ${task.cpus} ${meta.reference} ${fastq} > ${meta.sample_id}_reads2ref.sam
        samtools sort -o ${meta.sample_id}_reads2ref_sorted.bam ${meta.sample_id}_reads2ref.sam
		samtools index ${meta.sample_id}_reads2ref_sorted.bam
        """
}