process ALPHA_RAREFACTION {
    tag "Alpha rarefaction curve"
    publishDir "${params.outdir}/diversity/alpha-rarefaction", mode: 'copy'

    input:
    path table_qza
    path rooted_tree_qza

    output:
    path "alpha-rarefaction.qzv", emit: rarefaction_qzv

    script:
    """
    export NUMBA_CACHE_DIR=/tmp/numba_cache_\$\$
    export MPLCONFIGDIR=/tmp/mpl_cache_\$\$
    mkdir -p \$NUMBA_CACHE_DIR \$MPLCONFIGDIR

    qiime diversity alpha-rarefaction \
        --i-table ${table_qza} \
        --i-phylogeny ${rooted_tree_qza} \
        --p-max-depth ${params.max_rarefaction_depth} \
        --m-metadata-file ${params.metadata} \
        --o-visualization alpha-rarefaction.qzv
    """
}
