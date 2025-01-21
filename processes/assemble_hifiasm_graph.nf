process ASSEMBLE_HIFIASM {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/01_HIFIASM/", mode: "copy"

    input:
        tuple val(meta), path(fastq) // [ sample_id, fastq ]
    
    output:  
        tuple val(meta), path("${meta.sample_id}*.gfa"), emit: all_contigs
        tuple val(meta), path("${meta.sample_id}.bp.p_ctg.gfa"), emit: primary_contigs
    
    script:
        """
        hifiasm -t ${task.cpus} -o ${meta.sample_id} ${fastq}
        """	
}

