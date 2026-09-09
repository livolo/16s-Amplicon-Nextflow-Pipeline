process FASTP {

    tag "$sample_id"

    cpus 4

    container '/COLD_STORAGE/software/tools/fastp/fastp.sif'

    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(fwd), path(rev)

    output:
    tuple val(sample_id),
          path("${sample_id}_1.trimmed.fastq.gz"),
          path("${sample_id}_2.trimmed.fastq.gz"),
          emit: trimmed

    path "${sample_id}_fastp.json", emit: json

    path "${sample_id}_fastp.html", emit: html

    script:
    """
    fastp \
        -i ${fwd} \
        -I ${rev} \
        -o ${sample_id}_1.trimmed.fastq.gz \
        -O ${sample_id}_2.trimmed.fastq.gz \
        --json ${sample_id}_fastp.json \
        --html ${sample_id}_fastp.html \
        --thread ${task.cpus}
    """
}

