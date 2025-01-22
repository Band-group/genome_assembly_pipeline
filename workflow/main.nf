// Pipeline to assemble malaria genomes from long-read pacbio data

/* 
So far, the only input is a path to a samplesheet. There is currently no option to tweak with 
the settings of the processes (e.g. not able to change the hifiasm settings). Adding support for
this is a priority.
 */

// Processes
include { ASSEMBLE_HIFIASM                              } from '../processes/assemble_hifiasm_graph'
include { CREATE_FASTA                                  } from '../processes/create_fasta'
include { ALIGN_FA_TO_REF as ALIGN_FA_TO_REF_UNORIENTED } from '../processes/create_bam'
include { ALIGN_FA_TO_REF as ALIGN_FA_TO_REF_ORIENTED   } from '../processes/create_bam'
include { ORIENT_CONTIGS                                } from '../processes/orient_contigs'
include { RUN_QUAST_QC                                  } from '../processes/run_quast_qc'

// Subworkflows
include { PARSE_SAMPLESHEET                             } from '../subworkflows/parse_samplesheet'


workflow {
    reads_ch = PARSE_SAMPLESHEET ( params.samplesheet ) // pass params.input by reference. e.g. 
    reads_ch.map{it[0].sample_id}.view()

    ASSEMBLE_HIFIASM ( reads_ch )
    ch_primary_contigs = ASSEMBLE_HIFIASM.out.primary_contigs

    CREATE_FASTA ( ch_primary_contigs )
    ch_primary_fasta = CREATE_FASTA.out.contig_fa // => [ [sample_id, reference], fasta ]

    ALIGN_FA_TO_REF_UNORIENTED ( ch_primary_fasta )
    ch_primary_bam = ALIGN_FA_TO_REF_UNORIENTED.out.bam // => [ [sample_id, reference], bam ]

    // Join fasta and bam channels to ensure consistent order 
    ch_fasta_bam = ch_primary_fasta
        .join( ch_primary_bam, by: 0 ) // => [ meta, fasta, bam ]

    ORIENT_CONTIGS ( ch_fasta_bam )
    ch_oriented_ctgs = ORIENT_CONTIGS.out.oriented_fa
    
    ALIGN_FA_TO_REF_ORIENTED ( ch_oriented_ctgs )
    RUN_QUAST_QC ( ch_oriented_ctgs )
}

workflow.onComplete {
    // Create pipeline results directory
    file("${params.outdir}/results/pipeline_summary").mkdirs()

    // Write summary file
    def summaryFile = file("${params.outdir}/results/pipeline_summary/pipeline_summary.txt")
    summaryFile.text = """
        Pipeline finished successfully!
        Nextflow version: ${nextflow.version}
        Run time: ${workflow.duration}
        Run name: ${workflow.runName}
        Command line arguments: ${workflow.args.join(' ')}
        Parameters:
        ${params.inspect()}
        """
}