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
                [ [ sample_id:row.sample_id, reference:row.reference ] , row.fastq ]
            }
            .set { reads }

    emit:
        reads // channel: [ [meta], reads ]
}

/* 
Example CSV file:
sample_id,fastq,reference
s1,/path/to/file,/path/to/reference
s2,/path/to/file,/path/to/reference
s3,/path/to/file,/path/to/reference
s4,/path/to/file,/path/to/reference */