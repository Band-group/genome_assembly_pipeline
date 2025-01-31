process RUN_FASTQC {
    tag "$meta.sample_id"
    label "process_single"
    publishDir "${params.outdir}/results/01_FASTQ_QC_REPORTS/FastQC", mode: "copy"

    input:
        tuple val(meta), path(fastq) // [ sample_id, fastq ]

    output:
        path("*_fastqc.*"), emit: fastqc_report

    script:
        """
        fastqc ${fastq}
        """
}