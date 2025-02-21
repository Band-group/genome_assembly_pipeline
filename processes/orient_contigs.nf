def getSubstitutionString(meta){
	return [ // Each line is a sed string ('s/string_to_replace/new_string/')
		"'/Pf3D7_([0-9]+)_v3/${meta.sample_id}_\\1/'", // TO-DO: The strings to replace should not be hardcoded. They depend on the reference fa.
		"'/PF_apicoplast_genome_1/${meta.sample_id}_apicoplast/'",
		"'/Pf_M76611/${meta.sample_id}_mitchondrion/'"
	].join(' ')
}

process ORIENT_CONTIGS {
	tag "$meta.sample_id"
	label "process_single"
    publishDir "${params.outdir}/results/04_ORIENTED_CONTIGS/", mode: "copy"

	// requires iorek to be in path

    input:
		tuple val(meta), path(assembly_fasta), path(assembly_bam)

    output: 
        tuple val(meta), path("${meta.sample_id}_oriented.fa.gz"), emit: oriented_fa
		tuple val(meta), path("*.fai"), emit: oriented_fai

    script:
		def substitution = getSubstitutionString(meta)
        """
        order-orient \\
        	-alignments ${assembly_bam} \\
        	-fasta ${assembly_fasta} \\
        	-o ${meta.sample_id}_oriented.fa \\
			-substitute-contig-names ${substitution}

		bgzip ${meta.sample_id}_oriented.fa
		samtools faidx ${meta.sample_id}_oriented.fa.gz
        """
}