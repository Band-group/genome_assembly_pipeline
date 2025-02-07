// Subworkflow to subsample fastq files and create breadth of coverage plots

// Processes
include { GET_REF_LENGTH                              } from '../processes/get_ref_length'
include { SUBSAMPLE_FASTQ_READS                       } from '../processes/subsample_fastq_reads'
include { ALIGN_READS_TO_REF                          } from '../processes/align_reads_to_ref'
include { BAM_FILTER_ON_REGION                        } from '../processes/R/bam_filter_on_region'
include { COMPUTE_COVERAGE as COMPUTE_COVERAGE_REGION } from '../processes/compute_coverage'
include { COMPUTE_COVERAGE as COMPUTE_COVERAGE_ALL    } from '../processes/compute_coverage'
include { PLOT_COVERAGE                               } from '../processes/R/plot_coverage'

workflow SUBSAMPLE_FASTQ {
    take:
        ch_fastq         // => channel: [ [meta], reads ]

    main:
        GET_REF_LENGTH ( ch_fastq.map { it[0].reference } )
        ch_ref_length = GET_REF_LENGTH.out.ref_length

        SUBSAMPLE_FASTQ_READS ( ch_fastq, ch_ref_length )
        ch_subsampled_fastq = SUBSAMPLE_FASTQ_READS.out.sub_fastq
        
        ALIGN_READS_TO_REF ( ch_subsampled_fastq )
        if ( params.use_region ){ // Only run if use_regions is true
            BAM_FILTER_ON_REGION ( ALIGN_READS_TO_REF.out.aligned_reads_bam )
            ch_bams = BAM_FILTER_ON_REGION.out.reads_bam_filtered
        } else {
            ch_bams = ALIGN_READS_TO_REF.out.aligned_reads_bam
        }

        COMPUTE_COVERAGE ( ch_bams )
        ch_coverage_meta = COMPUTE_COVERAGE.out.coverage_data.collect{ list -> list[0] } // => [sample_id, reference, regions_bed]
        ch_coverage_data = COMPUTE_COVERAGE.out.coverage_data.collect{ list -> list[1] } // => [coverage_bed]
        
        PLOT_COVERAGE ( ch_coverage_meta, ch_coverage_data, ch_ref_length )
}