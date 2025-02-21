process RUN_FASTQC {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/01_FASTQ_QC_REPORTS/FastQC", mode: "copy"

    container "oras://community.wave.seqera.io/library/fastqc:0.12.1--104d26ddd9519960"

    input:
        tuple val(meta), path(fastq) // [ sample_id, fastq ]

    output:
        path("*_fastqc.*"), emit: fastqc_report

    script:
        """
        fastqc ${fastq}
        """
}