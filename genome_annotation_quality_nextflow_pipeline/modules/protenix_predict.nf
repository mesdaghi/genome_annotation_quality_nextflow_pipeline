#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process PROTENIX_PREDICT {
    tag { "${dataset_name}_${chunk_file.simpleName}" }

    cpus 8
    memory '16 GB'
    time '2d'
    clusterOptions '--gres=gpu:1'

    input:
    tuple val(dataset_name), path(chunk_file)

    output:
    tuple val(dataset_name), path("protenix_out/${dataset_name}"), emit: protenix_out

    script:
    """
    module load cuda/13.0.2
    module load python/3.10
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate protenix_env
    export PROTENIX_CACHE=/home/shahmes/protenix_cache

    mkdir -p protenix_out/${dataset_name}

    protenix predict \
        --input ${chunk_file} \
        --out_dir protenix_out/${dataset_name} \
        --seeds 101 \
        --model_name "protenix_mini_esm_v0.5.0" \
        --use_msa false \
        --sample 1
    """
}