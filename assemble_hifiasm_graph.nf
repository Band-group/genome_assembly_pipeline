process ASSEMBLE_HIFIASM {
    tag "$meta.sample_id"
    //label "process_single"
    publishDir "${projectDir}/results/01_HIFIASM/", mode: "copy"

    input:
        tuple val(meta), path(fastq) // [ sample_id, fastq ]
    
    output:
        //path "${sample_id}.txt"
        tuple val(meta), path("${meta.sample_id}*.gfa"), emit: all_contigs
        tuple val(meta), path("${meta.sample_id}.bp.p_ctg.gfa"), emit: primary_contigs
    
    script:
        """
		#echo ${sample_id} > ${sample_id}.txt
        hifiasm -t ${task.cpus} -o ${meta.sample_id} ${fastq}
        """	
}

