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