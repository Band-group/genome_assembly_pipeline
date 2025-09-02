process ALIGN_READS_TO_CONTIGS {
    tag "$meta.sample_id"
    label "process_low"
    publishDir "${params.outdir}/results/03_CONTIG_CORRECTION/reads2contigs/${meta.sample_id}", mode: "copy"

    container "oras://community.wave.seqera.io/library/minimap2_samtools:564e363a70ace2ad" // singularity

    input:
        tuple val(meta), path(fasta), path(fastq)

    output:
        tuple val(meta), path("${meta.sample_id}_reads2contig_sorted.bam"), emit: aligned_reads_bam
        tuple val(meta), path("*.bai"), emit: aligned_reads_bai
        file("versions.yml"), emit: versions

    script:
        """
        minimap2 -ax map-hifi -t ${task.cpus} ${fasta} ${fastq} > ${meta.sample_id}_reads2contig.sam
        samtools sort -o ${meta.sample_id}_reads2contig_sorted.bam ${meta.sample_id}_reads2contig.sam
		samtools index ${meta.sample_id}_reads2contig_sorted.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            minimap2: \$(minimap2 --version 2>&1)
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}