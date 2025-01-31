// Subworkflow to subsample fastq files and create breadth of coverage plots

// Processes
include { GET_REF_LENGTH           } from '../processes/get_ref_length'
include { SUBSAMPLE_FASTQ_READS    } from '../processes/subsample_fastq_reads'
include { ALIGN_READS_TO_REF       } from '../processes/align_reads_to_ref'
include { BAM_FILTER_ON_REGION     } from '../processes/R/bam_filter_on_region'
include { PLOT_BREADTH_OF_COVERAGE } from '../processes/plot_breadth_of_coverage'

workflow SUBSAMPLE_FASTQ {
    take:
        ch_fastq         // channel: [ [meta], reads ]
        regions_bed      // path to regions bed file

    main:
        Channel.fromPath( regions_bed )
            .set{ ch_regions_bed }

        ch_reference_length = GET_REF_LENGTH ( ch_fastq.map { it[0].reference } ).ref_length // get reference

        ch_subsampled_fastq = SUBSAMPLE_FASTQ_READS ( ch_fastq, ch_reference_length ).sub_fastq
        ch_aligned_reads = ALIGN_READS_TO_REF ( ch_subsampled_fastq ).aligned_reads_bam
        ch_region_restricted_bam = BAM_FILTER_ON_REGION ( ch_aligned_reads, ch_regions_bed ).reads_bam_filtered
        
        // Collect all BAMs before plotting
        PLOT_BREADTH_OF_COVERAGE ( ch_region_restricted_bam.collect(), ch_regions_bed )
}