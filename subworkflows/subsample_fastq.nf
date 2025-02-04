// Subworkflow to subsample fastq files and create breadth of coverage plots

// Processes
include { GET_REF_LENGTH                              } from '../processes/get_ref_length'
include { SUBSAMPLE_FASTQ_READS                       } from '../processes/subsample_fastq_reads'
include { ALIGN_READS_TO_REF                          } from '../processes/align_reads_to_ref'
include { BAM_FILTER_ON_REGION                        } from '../processes/R/bam_filter_on_region'
include { GET_CHROMOSOME_LENGTHS                      } from '../processes/get_chromosome_lengths'
include { COMPUTE_COVERAGE as COMPUTE_COVERAGE_REGION } from '../processes/compute_coverage'
include { COMPUTE_COVERAGE as COMPUTE_COVERAGE_ALL    } from '../processes/compute_coverage'

workflow SUBSAMPLE_FASTQ {
    take:
        ch_fastq         // => channel: [ [meta], reads ]

    main:
        GET_REF_LENGTH ( ch_fastq.map { it[0].reference } )
        ch_ref_length = GET_REF_LENGTH.out.ref_length

        SUBSAMPLE_FASTQ_READS ( ch_fastq, ch_ref_length )
        ch_subsampled_fastq = SUBSAMPLE_FASTQ_READS.out.sub_fastq
        
        ALIGN_READS_TO_REF ( ch_subsampled_fastq )
        ch_bams = ALIGN_READS_TO_REF.out.aligned_reads_bam

        if ( params.use_region ){ // Only run if use_regions is true
            BAM_FILTER_ON_REGION ( ch_bams )
            ch_region_filtered_bams = BAM_FILTER_ON_REGION.out.reads_bam_filtered
            GET_CHROMOSOME_LENGTHS ( ch_region_filtered_bams )
            ch_region_filtered_bams_cl = GET_CHROMOSOME_LENGTHS.out.chr_lengths // => [ [meta], bam, chr_lengths ]

            COMPUTE_COVERAGE_REGION ( ch_region_filtered_bams_cl )
        } else {
            // Add fake regions_bed 
            ch_bams 
                .map { row -> [ row[0], row[1], [] ] } 
                .set { ch_fregions_bams } // => [ [meta], bam, empty ]
            COMPUTE_COVERAGE_ALL ( ch_fregions_bams )
        }
}