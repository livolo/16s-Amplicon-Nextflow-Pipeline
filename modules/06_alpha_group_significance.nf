process ALPHA_GROUP_SIGNIFICANCE {
    tag "Alpha diversity group significance"
    publishDir "${params.outdir}/diversity/alpha-group-significance", mode: 'copy'

    input:
    path faith_pd_qza
    path shannon_qza
    path evenness_qza
    path observed_features_qza

    output:
    path "faith-pd-group-significance.qzv",         emit: faith_pd_sig
    path "shannon-group-significance.qzv",          emit: shannon_sig
    path "evenness-group-significance.qzv",         emit: evenness_sig
    path "observed-features-group-significance.qzv", emit: observed_features_sig

    script:
    """
    export NUMBA_CACHE_DIR=/tmp/numba_cache_\$\$
    export MPLCONFIGDIR=/tmp/mpl_cache_\$\$
    mkdir -p \$NUMBA_CACHE_DIR \$MPLCONFIGDIR

    qiime diversity alpha-group-significance \
        --i-alpha-diversity ${faith_pd_qza} \
        --m-metadata-file ${params.metadata} \
        --o-visualization faith-pd-group-significance.qzv

    qiime diversity alpha-group-significance \
        --i-alpha-diversity ${shannon_qza} \
        --m-metadata-file ${params.metadata} \
        --o-visualization shannon-group-significance.qzv

    qiime diversity alpha-group-significance \
        --i-alpha-diversity ${evenness_qza} \
        --m-metadata-file ${params.metadata} \
        --o-visualization evenness-group-significance.qzv

    qiime diversity alpha-group-significance \
        --i-alpha-diversity ${observed_features_qza} \
        --m-metadata-file ${params.metadata} \
        --o-visualization observed-features-group-significance.qzv
    """
}
