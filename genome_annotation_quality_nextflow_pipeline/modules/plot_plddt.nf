#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process PLOT_PLDDT {
    tag { dataset_name }


    publishDir "results/plots", mode: 'copy'

    input:
    tuple val(dataset_name), path(pkl_file)

    output:
    tuple val(dataset_name), path("plddt_*_${dataset_name}.png")

    script:
    """
    python ${projectDir}/bin/plot_plddt.py ${dataset_name}
    """
}