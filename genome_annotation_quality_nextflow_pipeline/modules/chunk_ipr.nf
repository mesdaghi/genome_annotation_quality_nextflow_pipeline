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