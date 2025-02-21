process ALIGN_READS_TO_GFA {
    tag "$meta.sample_id"
    label "process_high"
    publishDir "${params.outdir}/results/02_ASSEMBLIES/reads2gfa/${meta.sample_id}", mode: "copy"

    container "oras://community.wave.seqera.io/library/graphaligner:1.0.20--3fb1e479806451c9" // singularity

    input:
        tuple val(meta), path(assembly_graph), path(fastq) // => [ meta, primary contigs, fastq reads ]

    output:
        tuple val(meta), path("*gaf"), emit: reads2gfa_gaf
        tuple val(meta), path("*bam"), emit: reads2gfa_bam
    
    script: // TODO: Build a pangenome variant graph from all assemblies (using contigs or scaffolds? Probably makes most sense to use the contigs because scaffolds are reference-guided). Align reads from inidividual samples back to the graph.
        """
        GraphAligner -g ${assembly_graph} -f ${fastq} -t ${task.cpus} -a ${meta.sample_id}_reads2gfa.gaf
        """
}