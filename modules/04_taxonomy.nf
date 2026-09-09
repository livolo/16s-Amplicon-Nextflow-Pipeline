process TAXONOMY {
    tag "Taxonomic classification"
    publishDir "${params.outdir}/taxonomy", mode: 'copy'

    input:
    path rep_seqs_qza
    path classifier_qza
    path table_qza

    output:
    path "taxonomy.qza",       emit: taxonomy
    path "taxonomy.qzv",       emit: taxonomy_qzv
    path "taxa-bar-plots.qzv", emit: barplot

    script:
    """
    export NUMBA_CACHE_DIR=/tmp/numba_cache_\$\$
    export MPLCONFIGDIR=/tmp/mpl_cache_\$\$
    mkdir -p \$NUMBA_CACHE_DIR \$MPLCONFIGDIR

    qiime feature-classifier classify-sklearn \
        --i-classifier ${classifier_qza} \
        --i-reads ${rep_seqs_qza} \
        --p-n-jobs ${task.cpus} \
        --o-classification taxonomy.qza \
        --verbose

    qiime metadata tabulate \
        --m-input-file taxonomy.qza \
        --o-visualization taxonomy.qzv

    qiime taxa barplot \
        --i-table ${table_qza} \
        --i-taxonomy taxonomy.qza \
        --m-metadata-file ${params.metadata} \
        --o-visualization taxa-bar-plots.qzv
    """
}
