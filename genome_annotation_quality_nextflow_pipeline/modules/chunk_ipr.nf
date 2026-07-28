#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process CHUNK_IPR {
    tag { dataset_name }

    label 'cpu'
   
 
    input:
    tuple path(fasta_file), val(dataset_name)

    output:
    tuple val(dataset_name), path("ipr_chunks/*.fa")

    script:
    """
    mkdir -p ipr_chunks
    bash ${projectDir}/bin/balance_chunks.sh --input ${fasta_file} --chunks ${params.ipr_chunks} --outdir ipr_chunks
    """
}