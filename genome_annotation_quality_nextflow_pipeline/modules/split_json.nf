#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process SPLIT_JSON {
    cpus 1
    memory '2 GB'
    time '1h'
    tag { dataset_name }

    input:
    tuple path(json_file), val(dataset_name)

    output:
    tuple val(dataset_name), path("json_chunks/*.json")

    script:
    """
    mkdir -p json_chunks
    python ${projectDir}/bin/split_json.py ${json_file} ${params.chunk_size} json_chunks
    """
}