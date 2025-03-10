// Pipeline to assemble parasite genomes from long-read pacbio data.

// Processes
include { RUN_FASTQC                                    } from '../processes/run_fastqc'
include { ASSEMBLE_HIFIASM                              } from '../processes/assemble_hifiasm_graph'
include { ALIGN_READS_TO_GFA                            } from '../processes/align_reads_to_gfa'
include { ALIGN_READS_TO_REF                            } from '../processes/align_reads_to_ref'
include { CREATE_FA_FROM_GFA                            } from '../processes/create_fa_from_gfa'
include { ALIGN_FA_TO_REF as ALIGN_FA_TO_REF_UNORIENTED } from '../processes/align_fa_to_ref'
include { ALIGN_FA_TO_REF as ALIGN_FA_TO_REF_ORIENTED   } from '../processes/align_fa_to_ref'
include { ORIENT_CONTIGS                                } from '../processes/orient_contigs'
include { ALIGN_READS_TO_CONTIGS                        } from '../processes/align_reads_to_contigs'
include { COMPUTE_COVERAGE_PER_BASE                     } from '../processes/compute_coverage_per_base'
include { COMPUTE_COV_ZSCORES                           } from '../processes/compute_cov_zscores'
include { RUN_QUAST_QC                                  } from '../processes/run_quast_qc'

// Subworkflows
include { PARSE_SAMPLESHEET                             } from '../subworkflows/parse_samplesheet'
include { SUBSAMPLE_BAM                                 } from '../subworkflows/subsample_bam'
include { SCAFFOLD                                      } from '../subworkflows/scaffold'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// TODO: Use singularity so that individual processes run in dedicated containers.

workflow GENOMEASSEMBLER {
    ch_fastq = PARSE_SAMPLESHEET ( params.samplesheet )
    
    // QC
    RUN_FASTQC ( ch_fastq )
    SUBSAMPLE_BAM ( ch_fastq )
    
    // Assembly
    ASSEMBLE_HIFIASM ( ch_fastq )
    ch_primary_contigs = ASSEMBLE_HIFIASM.out.primary_contigs

    ch_graph_fastq = ch_primary_contigs 
        .join( ch_fastq, by: 0 ) // => [ [sample_id, reference], assembly_graph, fastq ]
    //ALIGN_READS_TO_GFA ( ch_graph_fastq )

    CREATE_FA_FROM_GFA ( ch_primary_contigs )
    ch_primary_fasta = CREATE_FA_FROM_GFA.out.contig_fa // => [ [sample_id, reference], contig_fa ]

    ALIGN_FA_TO_REF_UNORIENTED ( ch_primary_fasta )
    ch_primary_bam = ALIGN_FA_TO_REF_UNORIENTED.out.aligned_fa_bam // => [ [sample_id, reference], contig_bam ]

    // Join fasta and bam channels to ensure consistent order 
    ch_fasta_bam = ch_primary_fasta
        .join( ch_primary_bam, by: 0 ) // => [ meta, fasta, bam ]

    ORIENT_CONTIGS ( ch_fasta_bam )
    ch_oriented_ctgs = ORIENT_CONTIGS.out.oriented_fa // => [ meta, oriented_fasta ]
    ch_oriented_ctgs_meta_all = ORIENT_CONTIGS.out.oriented_fa.collect{ list -> list[0] } // => [sample_id, reference, regions_bed]
    ch_oriented_ctgs_data_all = ORIENT_CONTIGS.out.oriented_fa.collect{ list -> list[1] } // => [oriented_fasta]

    ch_oriented_ctgs_reads = ch_oriented_ctgs
        .join ( ch_fastq, by: 0 ) // => [ meta, oriented_fasta, fastq ]
    
    // TO-DO: below should be moved into a subworkflow
    // Align reads back to oriented contigs
    ALIGN_READS_TO_CONTIGS ( ch_oriented_ctgs_reads )
    ch_ctgs_reads_bam = ALIGN_READS_TO_CONTIGS.out.aligned_reads_bam
    // Compute per-base coverage and calculate window z-scores
    COMPUTE_COVERAGE_PER_BASE ( ch_ctgs_reads_bam )
    ch_perbase_cov = COMPUTE_COVERAGE_PER_BASE.out.coverage_data
    COMPUTE_COV_ZSCORES ( ch_perbase_cov )
    
    ALIGN_FA_TO_REF_ORIENTED ( ch_oriented_ctgs )
    RUN_QUAST_QC ( ch_oriented_ctgs )
}