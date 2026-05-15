#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process PLOT_PLDDT {
    tag { dataset_name }

    label 'cpu'

    publishDir "results/plots", mode: 'copy' // having dirs of the same name sometimes causes issues with singularity.

    input:
    tuple val(dataset_name), path(pkl_file)

    output:
    tuple val(dataset_name), path("plddt_*_${dataset_name}.png")

    script:
    """
    export MPLCONFIGDIR=\$(mktemp -d)
    
    python3 ${projectDir}/bin/plot_plddt.py ${dataset_name}
    """
}