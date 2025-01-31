process ALIGN_READS_TO_GFA {
    tag "$meta.sample_id"
    label "process_low"
    publishDir "${params.outdir}/results/02_READS_ALIGNED_TO_GFA/", mode: "copy"

    input:
        tuple val(meta), path(assembly_graph)
        tuple val(meta2), path(fastq)

    output:
        tuple val(meta), path("*gam"), emit: reads2gfa_gam
        tuple val(meta), path("*bam"), emit: reads2gfa_bam
    
    script:
        """
        GraphAligner -g ${assembly_graph} -f ${fastq} -t ${task.cpus} -a ${meta.sample_id}_reads2gfa.gam
        """


    
}