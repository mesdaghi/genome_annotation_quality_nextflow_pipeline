nextflow.enable.dsl=2

// ==================== Processes ==================== //

include { FASTA_TO_JSON } from './modules/fasta_to_json.nf'
include {SPLIT_JSON} from './modules/split_json.nf'
include { PROTENIX_PREDICT } from './modules/protenix_predict.nf'
include { COLLECT_CHUNKS } from './modules/collect_chunks.nf'
include { PROCESS_MODELS } from './modules/process_models.nf'
include { PLOT_PLDDT } from './modules/plot_plddt.nf'
include { METAPREDICT_DISORDER } from './modules/metapredict_disorder.nf'
include { PLOT_METAPREDICT } from './modules/plot_metapredict.nf'
include { PSAURON_RUN } from './modules/psauron_run.nf'
include { PLOT_PSAURON } from './modules/plot_psauron.nf'
include { GENERATE_REPORT } from './modules/generate_report.nf'

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

    // ================= REPORT STAGE ================= //
    // Collect all PNGs per dataset_name into a single list, then generate report

    report_ch = plddt_plot_ch
        .join(metapredict_plot_ch, by: 0)
        .join(psauron_plot_ch, by: 0)
        .map { dataset_name, plddt_png, metapredict_pngs, psauron_pngs ->
            // Flatten into a single list regardless of how many PNGs each process emits
            def pngs = []
            pngs += (plddt_png instanceof List ? plddt_png : [plddt_png])
            pngs += (metapredict_pngs instanceof List ? metapredict_pngs : [metapredict_pngs])
            pngs += (psauron_pngs instanceof List ? psauron_pngs : [psauron_pngs])
            tuple(dataset_name, pngs)
        }

    GENERATE_REPORT(report_ch)

}
