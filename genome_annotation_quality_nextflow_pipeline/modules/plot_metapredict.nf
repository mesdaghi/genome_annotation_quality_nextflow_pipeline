#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process PLOT_METAPREDICT {
    tag { dataset_name }


    publishDir "results/metapredict_plots", mode: 'copy'

    input:
    tuple val(dataset_name), path(csv_file)

    output:
    tuple val(dataset_name), path("*.png")

    script:
    """
    python ${projectDir}/bin/plot_metapredict.py ${csv_file}

    # Where is the output file goign.
    """
}