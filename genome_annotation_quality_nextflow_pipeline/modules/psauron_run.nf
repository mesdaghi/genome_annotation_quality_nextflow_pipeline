#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process PSAURON_RUN {
    tag { dataset_name }

   
    publishDir "results/psauron", mode: 'copy'
    clusterOptions '--gres=gpu:1'

    input:
    tuple path(fasta_file), val(dataset_name)

    output:
    tuple val(dataset_name), path("${dataset_name}_psauron.csv")

    script:
    """
    echo "Running on: \$(hostname)"
    echo "Start time: \$(date)"

    module load cuda/13.0.2

    nvidia-smi

    psauron \
        -i ${fasta_file} \
        -o ${dataset_name}_psauron.csv \
        -p

    echo "End time: \$(date)"
    """
}