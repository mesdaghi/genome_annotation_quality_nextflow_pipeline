#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process PLOT_PSAURON {
    tag { dataset_name }

 

    publishDir "results/psauron_plots", mode: 'copy' // having dirs of the same name sometimes causes issues with singularity.

    input:
    tuple val(dataset_name), path(psauron_csv)

    output:
    tuple val(dataset_name), path("*.png")

    script:
    """
    python3 ${projectDir}/bin/plot_psauron_distribution.py \
        ${projectDir}/bin/combined_psauron_results.csv \ # This shoud not be hardcoded...
        ${psauron_csv}
    """
}