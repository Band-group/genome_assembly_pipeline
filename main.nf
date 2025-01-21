// Pipeline to assemble malaria genomes from long-read pacbio data

/* 
So far, the only input is a path to a samplesheet. There is currently no option to tweak with 
the settings of the processes (e.g. not able to change the hifiasm settings). Adding support for
this is a priority.
 */

// Processes
include { ASSEMBLE_HIFIASM                 } from './assemble_hifiasm_graph'
include { CREATE_FASTA                     } from './create_fasta'
include { CREATE_BAM                       } from './create_bam'
include { CREATE_BAM as CREATE_BAM_CONTIGS } from './create_bam'
include { ORIENT_CONTIGS                   } from './orient_contigs'

// Subworkflows
include { PARSE_SAMPLESHEET                } from './parse_samplesheet'


workflow {
    //reads_ch = PARSE_SAMPLESHEET ( params.samplesheet ) // pass params.input by reference. e.g. 
    //reads_ch.map{it[0].sample_id}.view()

    //ASSEMBLE_HIFIASM ( reads_ch )
    //ch_primary_contigs = ASSEMBLE_HIFIASM.out.primary_contigs
    
    //CREATE_FASTA ( ch_primary_contigs )
    ch_primary_contigs = Channel
        .from( params.primary_contigs )
        .map { row ->
            tuple( row[0], file(row[1]) )
        }
    CREATE_FASTA ( ch_primary_contigs )
    ch_primary_fasta = CREATE_FASTA.out.contig_fa // => [ [sample_id, reference], fasta ]

    CREATE_BAM ( ch_primary_fasta )
    ch_primary_bam = CREATE_BAM.out.bam // => [ [sample_id, reference], bam ]

    // Join fasta and bam channels to ensure consistent order 
    ch_fasta_bam = ch_primary_fasta
        .join( ch_primary_bam, by: 0 ) // => [ meta, fasta, bam ]

    ORIENT_CONTIGS ( ch_fasta_bam )
    ch_oriented_ctgs = ORIENT_CONTIGS.out.oriented_fa
    
    CREATE_BAM_CONTIGS ( ch_oriented_ctgs )
}