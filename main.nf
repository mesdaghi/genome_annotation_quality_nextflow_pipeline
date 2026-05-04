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
    tuple val(dataset_name), path("plddt_*_${dataset_name}.png")

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
    tuple val(dataset_name), path("*.png")

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
    tuple val(dataset_name), path("*.png")

    script:
    """
    python ${projectDir}/bin/plot_psauron_distribution.py \
        ${projectDir}/bin/combined_psauron_results.csv \
        ${psauron_csv}
    """
}

process GENERATE_REPORT {
    tag { dataset_name }

    cpus 1
    memory '2 GB'
    time '30m'

    publishDir "results", mode: 'copy'

    input:
    tuple val(dataset_name), path(all_pngs)

    output:
    path "${dataset_name}_report.pdf"

    script:
    """
    pip install reportlab pillow --quiet

    python ${projectDir}/bin/generate_report.py \
        ${dataset_name} \
        ${dataset_name}_report.pdf \
        ${all_pngs}
    """
}


// ==================== InterPro Branch ==================== //

process CHUNK_IPR {
    tag { dataset_name }

    cpus 1
    memory '2 GB'
    time '1h'

    input:
    tuple path(fasta_file), val(dataset_name)

    output:
    tuple val(dataset_name), path("ipr_chunks/*.fasta")

    script:
    """
    mkdir -p ipr_chunks
    bash ${projectDir}/bin/chunk_ipr.sh ${fasta_file} 20 ipr_chunks
    """
}

process INTERPROSCAN {
    tag { "${dataset_name}_${chunk_file.simpleName}" }

    cpus 16
    memory '64 GB'
    time '24h'

    input:
    tuple val(dataset_name), path(chunk_file)

    output:
    tuple val(dataset_name), path("*.tsv"), path("*.xml")

    script:
    """
    module purge
    module load adoptopenjdk/11.0.12+7
    module load interproscan/5.66-98.0

    interproscan.sh \
        -i ${chunk_file} \
        -f TSV,XML \
        --cpu ${task.cpus}
    """
}

process CONCAT_IPR {
    tag { dataset_name }

    cpus 2
    memory '8 GB'
    time '2h'

    publishDir "results/interpro", mode: 'copy'

    input:
    tuple val(dataset_name), path(tsv_files), path(xml_files)

    output:
    tuple val(dataset_name), path("concatenated.tsv"), path("merged_output.xml")

    script:
    """
    python ${projectDir}/bin/concat_xml.py ${xml_files} -o merged_output.xml
    python ${projectDir}/bin/concat_tsv.py ${tsv_files}
    """
}

process PLOT_INTERPRO {
    tag { dataset_name }

    cpus 2
    memory '4 GB'
    time '2h'

    publishDir "results/interpro_plots", mode: 'copy'

    input:
    tuple val(dataset_name), path(fasta_file), path(xml_file)

    output:
    tuple val(dataset_name),
          path("interpro_domain_coverage.png"),
          path("interpro_merged_coverage_distribution.png"),
          path("interpro_summary.tsv")

    script:
    """
    python ${projectDir}/bin/plot_interpro.py \
        ${fasta_file} \
        ${xml_file} \
        --json ${params.ipr_reference_json} \
        --name "${dataset_name}" \
        --outdir .

    # Rename to stable filenames so the report step can find them
    mv query_ipr_coverage.png                       interpro_domain_coverage.png
    mv query_ipr_merged_coverage_distribution.png   interpro_merged_coverage_distribution.png
    mv query_ipr_summary.tsv                        interpro_summary.tsv
    """
}

// params.run_interpro behaviour:
//   false (default)        → interpro branch disabled entirely
//   true (bare flag)       → full branch: chunk → interproscan → concat → plot
//   "/path/to/file.xml"    → plot only, using the provided InterProScan XML
//                            (paired with the dataset's FASTA from fasta_ch)

params.run_interpro = false

// Path to the precomputed reference dataset of model organisms used for
// the comparison plots. See bin/ipr_coverage.README.md for how to build it.
params.ipr_reference_json = "${projectDir}/bin/ipr_coverage.json"

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

    plddt_plot_ch = PLOT_PLDDT(pkl_ch)

    // ================= METAPREDICT BRANCH ================= //

    metapredict_ch = METAPREDICT_DISORDER(fasta_ch)

    metapredict_plot_ch = PLOT_METAPREDICT(metapredict_ch)

    // ================= PSAURON BRANCH ================= //

    psauron_ch = PSAURON_RUN(fasta_ch)

    psauron_plot_ch = PLOT_PSAURON(psauron_ch)

    // ================= INTERPRO BRANCH (optional) ================= //

    if ( params.run_interpro instanceof String ) {

        // User supplied a pre-existing XML — plot only, paired with the
        // dataset's own FASTA so the new script can compute totals.
        def supplied_xml = file(params.run_interpro, checkIfExists: true)
        if ( supplied_xml.extension.toLowerCase() != 'xml' ) {
            error "params.run_interpro must point to an InterProScan XML " +
                  "file (got: ${supplied_xml})"
        }
        interpro_plot_ch = fasta_ch
            .map { fasta_file, dataset_name ->
                tuple(dataset_name, fasta_file, supplied_xml)
            }
        interpro_plot_ch = PLOT_INTERPRO(interpro_plot_ch)

    } else if ( params.run_interpro == true ) {

        // Full branch
        ipr_chunks_ch = CHUNK_IPR(fasta_ch)
            .flatMap { dataset_name, chunks ->
                chunks.collect { chunk_file -> tuple(dataset_name, chunk_file) }
            }

        ipr_results_ch = INTERPROSCAN(ipr_chunks_ch)

        ipr_collected_ch = ipr_results_ch
            .groupTuple(by: 0)
            .map { dataset_name, tsv_list, xml_list ->
                def tsvs = tsv_list instanceof List ? tsv_list.flatten() : [tsv_list]
                def xmls  = xml_list  instanceof List ? xml_list.flatten()  : [xml_list]
                tuple(dataset_name, tsvs, xmls)
            }

        ipr_concat_ch = CONCAT_IPR(ipr_collected_ch)

        // Pair the merged XML back up with each dataset's FASTA for plotting
        plot_input_ch = ipr_concat_ch
            .map { dataset_name, tsv, xml -> tuple(dataset_name, xml) }
            .join(
                fasta_ch.map { fasta_file, dataset_name -> tuple(dataset_name, fasta_file) },
                by: 0
            )
            .map { dataset_name, xml, fasta_file -> tuple(dataset_name, fasta_file, xml) }

        interpro_plot_ch = PLOT_INTERPRO(plot_input_ch)

    } else {

        // Disabled — empty channel so the report join is unaffected
        interpro_plot_ch = Channel.empty()

    }

    // ================= REPORT STAGE ================= //
    // Collect all PNGs per dataset_name into a single list, then generate report

    report_ch = plddt_plot_ch
        .join(metapredict_plot_ch, by: 0)
        .join(psauron_plot_ch, by: 0)

    if ( params.run_interpro != false ) {
        report_ch = report_ch
            .join(interpro_plot_ch, by: 0)
            .map { dataset_name, plddt_png, metapredict_pngs, psauron_pngs,
                   interpro_bar_png, interpro_dist_png, interpro_tsv ->
                def pngs = []
                pngs += (plddt_png instanceof List ? plddt_png : [plddt_png])
                pngs += (metapredict_pngs instanceof List ? metapredict_pngs : [metapredict_pngs])
                pngs += (psauron_pngs instanceof List ? psauron_pngs : [psauron_pngs])
                // Both interpro PNGs go into the report; the TSV is published
                // to results/interpro_plots but not embedded.
                pngs += [interpro_bar_png, interpro_dist_png]
                tuple(dataset_name, pngs)
            }
    } else {
        report_ch = report_ch
            .map { dataset_name, plddt_png, metapredict_pngs, psauron_pngs ->
                def pngs = []
                pngs += (plddt_png instanceof List ? plddt_png : [plddt_png])
                pngs += (metapredict_pngs instanceof List ? metapredict_pngs : [metapredict_pngs])
                pngs += (psauron_pngs instanceof List ? psauron_pngs : [psauron_pngs])
                tuple(dataset_name, pngs)
            }
    }

    GENERATE_REPORT(report_ch)

}



