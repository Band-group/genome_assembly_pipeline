process ASSEMBLE_HIFIASM {
    tag "$meta.sample_id"
    label "process_medium"
    publishDir "${params.outdir}/results/02_ASSEMBLIES/hifiasm/${meta.sample_id}_nhaps${params.hifiasm_n_hap}_purgelvl${params.hifiasm_purge_dup_lvl}_D${params.hifiasm_d}_N${params.hifiasm_n}", mode: "copy"

    container "oras://community.wave.seqera.io/library/hifiasm:0.24.0--8d78d9b82cf20802" // singularity

    input:
        tuple val(meta), path(fastq) // [ sample_id, fastq ]
    
    output:  
        tuple val(meta), path("${meta.sample_id}.bp.p_ctg.gfa"), emit: primary_contigs
        tuple val(meta), path("${meta.sample_id}*"), emit: all_output
        path("versions.yml"), emit: versions
    
    script:
        def n_hap = params.hifiasm_n_hap ? "--n-hap ${params.hifiasm_n_hap}" : ''
        def purge_dup_lvl = "-l${params.hifiasm_purge_dup_lvl}"
        def d = "-D ${params.hifiasm_d}"
        def n = "-N ${params.hifiasm_n}"

        """
        hifiasm -t ${task.cpus} -o ${meta.sample_id} ${n_hap} ${purge_dup_lvl} ${d} ${n} ${fastq}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hifiasm: \$(hifiasm --version 2>&1)
        END_VERSIONS
        """	
}

