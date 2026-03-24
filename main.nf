nextflow.enable.dsl=2

// ==================== Processes ==================== //

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

process COLLECT_CHUNKS {
    tag { dataset_name }

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

process PROCESS_MODELS {
    tag { dataset_name }

    cpus 4
    memory '8 GB'
    time '4h'

    input:
    tuple val(dataset_name), path(pred_dir)

    output:
    tuple val(dataset_name), path("plddt_all_values_${dataset_name}_all_one.pkl")

    script:
    """
    python ${projectDir}/bin/process_models.py \
      ${pred_dir} \
      plddt_all_values_${dataset_name}_all_one.pkl
    """
}

process PLOT_PLDDT {
    tag { dataset_name }

    cpus 2
    memory '4 GB'
    time '2h'
    publishDir "results/plots", mode: 'copy'

    input:
    tuple val(dataset_name), path(pkl_file)

    output:
    path "*.png"

    script:
    """
    python ${projectDir}/bin/plot_plddt.py ${dataset_name}
    """
}

process METAPREDICT_DISORDER {
    tag { dataset_name }

    cpus 8
    memory '16 GB'
    time '2d'

    publishDir "results/metapredict", mode: 'copy'

    // GPU request (SLURM-style)
    clusterOptions '--gres=gpu:1'

    input:
    tuple path(fasta_file), val(dataset_name)

    output:
    tuple val(dataset_name), path("${dataset_name}_metapredict.csv")

    script:
    """
    echo "Running on: \$(hostname)"
    echo "Start time: \$(date)"

    module load cuda/13.0.2
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate metapredict

    nvidia-smi

    metapredict-predict-disorder \
        -o ${dataset_name}_metapredict.csv \
        ${fasta_file}

    echo "End time: \$(date)"
    """
}

process PLOT_METAPREDICT {
    tag { dataset_name }

    cpus 2
    memory '4 GB'
    time '2h'

    publishDir "results/metapredict_plots", mode: 'copy'

    input:
    tuple val(dataset_name), path(csv_file)

    output:
    path "*.png"

    script:
    """
    python ${projectDir}/bin/plot_metapredict.py ${csv_file}
    """
}

process PSAURON_RUN {
    tag { dataset_name }

    cpus 8
    memory '16 GB'
    time '2d'

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

process PLOT_PSAURON {
    tag { dataset_name }

    cpus 2
    memory '4 GB'
    time '2h'

    publishDir "results/psauron_plots", mode: 'copy'

    input:
    tuple val(dataset_name), path(psauron_csv)

    output:
    path "*.png"

    script:
    """
    python ${projectDir}/bin/plot_psauron_distribution.py \
        ${projectDir}/bin/combined_psauron_results.csv \
        ${psauron_csv}
    """
}



workflow {

    fasta_ch = Channel
        .fromPath(params.fasta, checkIfExists: true)
        .map { fasta -> tuple(fasta, fasta.baseName) }

    // ================= PROTENIX BRANCH ================= //

    json_ch = FASTA_TO_JSON(fasta_ch)

    split_ch = SPLIT_JSON(json_ch)

    chunks_ch = split_ch
        .flatMap { dataset_name, chunks ->
            chunks.collect { chunk_file -> tuple(dataset_name, chunk_file) }
        }

    pred_ch = PROTENIX_PREDICT(chunks_ch)

    collected_ch = COLLECT_CHUNKS(pred_ch)

    pkl_ch = PROCESS_MODELS(collected_ch)

    PLOT_PLDDT(pkl_ch)

    // ================= METAPREDICT BRANCH ================= //

    metapredict_ch = METAPREDICT_DISORDER(fasta_ch)

    PLOT_METAPREDICT(metapredict_ch)


    // ================= PSAURON BRANCH ================= //

    psauron_ch = PSAURON_RUN(fasta_ch)

    PLOT_PSAURON(psauron_ch)

}

