#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process PLOT_PSAURON {
    tag { dataset_name }

    cpus 2
    memory '4 GB'
    time '2h'

    publishDir "results/psauron_plots", mode: 'copy'

    input:
    tuple val(dataset_name), path(psauron_csv)

    output:
    tuple val(dataset_name), path("*.png")

    script:
    """
    python ${projectDir}/bin/plot_psauron_distribution.py \
        ${projectDir}/bin/combined_psauron_results.csv \ # This shoud not be hardcoded...
        ${psauron_csv}
    """
}