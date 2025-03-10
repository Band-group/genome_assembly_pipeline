// Subsample a bam file to a target depth, discarding secondary alignments.

process SUBSAMPLE_BAM_READS {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/01_FASTQ_QC_REPORTS/coverage/subsampled_bams/${meta.sample_id}", mode: "copy"

    container "oras://community.wave.seqera.io/library/samtools:1.21--84c9d77c3901e90b" // singularity

    input:
        tuple val(meta), path(aligned_reads)
        val(reference_length)

    output:
        tuple val(meta), path("${meta.sample_id}_targetdepth_${params.target_depth}_subsampled.bam"), emit: sub_bam

    script:
        def target_depth = "${params.target_depth}"
        """
        TARGET_BASES=\$(( ${reference_length} * ${target_depth} ))
        
        {
            samtools view -H ${aligned_reads};
            samtools view -F 260 ${aligned_reads} | shuf;
        } | awk -v target=\$TARGET_BASES '
            BEGIN { total = 0 }
            /^@/ { print; next } # print header lines
            {
                rec = \$0
                len = length(\$10)
                total += len
                print rec

                if(total >= target) { exit }
            }
        ' | samtools view -b > ${meta.sample_id}_targetdepth_${target_depth}_subsampled.bam
        """
}