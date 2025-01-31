process SUBSAMPLE_FASTQ_READS {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/01_FASTQ_QC_REPORTS/coverage/subsampled_fastq", mode: "copy"

    input:
        tuple val(meta), path(fastq)
        val(reference_length)

    output:
        tuple val(meta), path("${meta.sample_id}_subsampled.fastq.gz"), emit: sub_fastq
    
    script:
        def target_depth = "${params.target_depth}"
        """
        # Calculate total bases in fastq
        BASES=\$(zcat ${fastq} | paste - - - - | cut -f2 | tr -d '\n' | wc -c)
    
        # Calculate current and target coverage
        CURRENT_COVERAGE=\$(( \$BASES / ${reference_length} ))
        TARGET_COVERAGE=\$(( ${reference_length} * ${params.target_depth} ))
        
        # Calculate sampling fraction
        FRACTION=\$(awk "BEGIN {f=${params.target_depth}/\$CURRENT_COVERAGE; print (f>1?1:f)}")
        
        # Subsample fastq file
        seqtk sample -s100 ${fastq} \$FRACTION | gzip > ${meta.sample_id}_\${FRACTION}_subsampled.fastq.gz
        """
}