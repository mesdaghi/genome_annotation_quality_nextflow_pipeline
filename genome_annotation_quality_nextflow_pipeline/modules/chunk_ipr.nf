#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process CHUNK_IPR {
    tag { dataset_name }

 
    input:
    tuple path(fasta_file), val(dataset_name)

    output:
    tuple val(dataset_name), path("ipr_chunks/*.fasta")

    script:
    """
    mkdir -p ipr_chunks
    bash ${projectDir}/bin/chunk_ipr.sh ${fasta_file} 20 ipr_chunks
    """
}