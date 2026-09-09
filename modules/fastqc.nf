process FASTQC {

    tag "$sample_id"

    cpus 4

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"),  emit: zip

    script:
    """
    fastqc \
        -t ${task.cpus} \
           ${reads}
    """
}
