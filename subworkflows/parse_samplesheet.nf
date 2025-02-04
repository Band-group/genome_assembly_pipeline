// Subworkflow to output a channel

workflow PARSE_SAMPLESHEET {
    take:
        samplesheet

    main:
        Channel.fromPath(samplesheet)
            .splitCsv( header: true, sep:"," )
            .map { row ->
                /* if ( row.type != 'pacbio' && row.type != 'nanopore' ) {
                    exit 1, "ERROR: ${row.type} is not a valid type. Valid types are \'nanopore\' or \'pacbio\'."
                } */
                if ( !file(row.fastq).exists() ) {
                    exit 1, "ERROR: FastQ file not found: ${row.fastq}"
                }
                if ( !file(row.reference).exists() ){
                    exit 1, "ERROR: Reference assembly not found: ${row.reference}"
                }
                // If 'regions_bed' is empty, set it to null
                if ( row.regions_bed.trim() == '' ) {
                    row.regions_bed = null
                } else if ( !file(row.regions_bed).exists() ){
                    exit 1, "ERROR: BED file not found: ${row.regions_bed}"
                }
                [ [ sample_id:row.sample_id, reference:row.reference, regions_bed:row.regions_bed ] , row.fastq ]
            }
            .set { reads }

    emit:
        reads // channel: [ [meta], reads ]
}

/* 
Example CSV file:
sample_id,fastq,reference,regions_bed
s1,/path/to/file,/path/to/reference,path/to/regions_bed
s2,/path/to/file,/path/to/reference,path/to/regions_bed
s3,/path/to/file,/path/to/reference,path/to/regions_bed
s4,/path/to/file,/path/to/reference,path/to/regions_bed 

regions_bed is not required. If you don't want to calculate coverage statistics for specific regions then either leave the column empty or set it to null.
*/