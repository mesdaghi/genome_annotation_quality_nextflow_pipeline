#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process COLLECT_CHUNKS {
    tag { dataset_name }

    label 'cpu'
   
    input:
    tuple val(dataset_name), path(pred_dir)

    output:
    tuple val(dataset_name), path("${dataset_name}_all_predictions")

    script:
    """
    mkdir -p ${dataset_name}_all_predictions
    cp -r ${pred_dir}/* ${dataset_name}_all_predictions/
    """
}