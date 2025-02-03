#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENOMEASSEMBLER } from "./workflow/genomeassembler"

workflow {
    GENOMEASSEMBLER ()
}
