process ASSEMBLE_HIFIASM {
    tag "$meta.sample_id"
    label "process_medium"
    publishDir "${params.outdir}/results/01_HIFIASM/", mode: "copy"

    input:
        tuple val(meta), path(fastq) // [ sample_id, fastq ]
    
    output:  
        tuple val(meta), path("${meta.sample_id}*.gfa"), emit: all_contigs
        tuple val(meta), path("${meta.sample_id}.bp.p_ctg.gfa"), emit: primary_contigs
    
    script:
        def n_haps = params.hifiasm_n_haps ? "--n-haps ${params.hifiasm_n_haps}" : '' // 1 for haploid assemblies
        def purge_dup_lvl = params.hifiasm_purge_dup_lvl ? "-l ${params.hifiasm_purge_dup_lvl}" : ''

        """
        hifiasm -t ${task.cpus} -o ${meta.sample_id} $n_haps $purge_dup_lvl ${fastq}
        """	
}

