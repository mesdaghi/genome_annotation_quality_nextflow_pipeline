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
include { CHUNK_IPR } from './modules/chunk_ipr.nf'
include { INTERPROSCAN } from './modules/interproscan.nf'
include { CONCAT_IPR } from './modules/concat_ipr.nf'
include { PLOT_INTERPRO } from './modules/plot_interpro.nf'


workflow {

    fasta_ch = Channel
        .fromPath(params.fasta, checkIfExists: true)
        .map { fasta -> tuple(fasta, fasta.baseName) }

    // ================= PROTENIX BRANCH ================= //

    json_ch = FASTA_TO_JSON(fasta_ch)

    split_ch = SPLIT_JSON(json_ch)
    chunks_ch = split_ch.transpose()


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
        ipr_chunks_ch = CHUNK_IPR(fasta_ch).transpose()

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
