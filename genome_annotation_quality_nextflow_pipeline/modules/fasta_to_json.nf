#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process FASTA_TO_JSON {
    cpus 1
    memory '2 GB'
    time '1h'
    tag { dataset_name }

    input:
    tuple path(fasta_file), val(dataset_name)

    output:
    tuple path("out.json"), val(dataset_name)

    script:
    """
    python ${projectDir}/bin/fasta_to_protenix_json_all.py ${fasta_file} out.json
    """
}