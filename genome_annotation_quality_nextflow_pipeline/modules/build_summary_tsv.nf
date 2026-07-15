#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * BUILD_SUMMARY_TSV — writes one TSV per dataset with the five headline
 * metrics shown on the pipeline's plots.
 *
 * Two process variants because Nextflow's `path` input doesn't have a
 * clean "optional" mode in DSL2: WITH_IPR takes the FASTA + XML for the
 * InterPro hit-rate column; the base variant leaves it blank. main.nf
 * picks the right one based on params.run_interpro.
 */


process BUILD_SUMMARY_TSV {
    tag { dataset_name }

    label 'cpu'

    publishDir "results", mode: 'copy'

    input:
    tuple val(dataset_name),
          path(plddt_pkl),
          path(metapredict_csv),
          path(psauron_csv)

    output:
    tuple val(dataset_name), path("${dataset_name}_summary.tsv")

    script:
    """
    export MPLCONFIGDIR=\$(mktemp -d)

    python3 ${projectDir}/bin/build_summary_tsv.py \\
        --name "${dataset_name}" \\
        --plddt-pkl ${plddt_pkl} \\
        --plddt-reference-csv ${projectDir}/reference/plddt_model_organisms.csv \\
        --metapredict-csv ${metapredict_csv} \\
        --psauron-csv ${psauron_csv} \\
        --out ${dataset_name}_summary.tsv
    """
}


process BUILD_SUMMARY_TSV_WITH_IPR {
    tag { dataset_name }

    label 'cpu'

    publishDir "results", mode: 'copy'

    input:
    tuple val(dataset_name),
          path(plddt_pkl),
          path(metapredict_csv),
          path(psauron_csv),
          path(ipr_fasta),
          path(ipr_xml)

    output:
    tuple val(dataset_name), path("${dataset_name}_summary.tsv")

    script:
    """
    export MPLCONFIGDIR=\$(mktemp -d)

    python3 ${projectDir}/bin/build_summary_tsv.py \\
        --name "${dataset_name}" \\
        --plddt-pkl ${plddt_pkl} \\
        --plddt-reference-csv ${projectDir}/reference/plddt_model_organisms.csv \\
        --metapredict-csv ${metapredict_csv} \\
        --psauron-csv ${psauron_csv} \\
        --ipr-fasta ${ipr_fasta} \\
        --ipr-xml ${ipr_xml} \\
        --plot-interpro-dir ${projectDir}/bin \\
        --out ${dataset_name}_summary.tsv
    """
}

