process DIFFERENTIAL_ABUNDANCE {
    tag "ANCOM-BC differential abundance"
    publishDir "${params.outdir}/differential-abundance", mode: 'copy'

    input:
    path table_qza
    path taxonomy_qza

    output:
    path "ancombc-health-state.qza",                        emit: ancombc_results
    path "ancombc-health-state-barplot.qzv",                emit: ancombc_barplot
    path "ancombc-health-state-barplot-relaxed.qzv",        emit: ancombc_barplot_relaxed
    path "exported-asv",                                    emit: ancombc_stats_dir
    path "table-l6.qza",                                    emit: table_genus
    path "ancombc-genus-health-state.qza",                  emit: ancombc_genus_results
    path "ancombc-genus-health-state-barplot.qzv",          emit: ancombc_genus_barplot
    path "ancombc-genus-health-state-barplot-relaxed.qzv",  emit: ancombc_genus_barplot_relaxed
    path "exported-genus",                                  emit: ancombc_genus_stats_dir

    script:
    """
    export NUMBA_CACHE_DIR=/tmp/numba_cache_\$\$
    export MPLCONFIGDIR=/tmp/mpl_cache_\$\$
    mkdir -p \$NUMBA_CACHE_DIR \$MPLCONFIGDIR

    # ASV-level differential abundance
    qiime composition ancombc \
        --i-table ${table_qza} \
        --m-metadata-file ${params.metadata} \
        --p-formula 'Health_State' \
        --o-differentials ancombc-health-state.qza \
        --verbose

    # Strict barplot (standard threshold, for formal reporting)
    qiime composition da-barplot \
        --i-data ancombc-health-state.qza \
        --p-significance-threshold 0.05 \
        --o-visualization ancombc-health-state-barplot.qzv

    # Relaxed barplot (exploratory only - NOT for formal significance claims)
    qiime composition da-barplot \
        --i-data ancombc-health-state.qza \
        --p-significance-threshold 0.2 \
        --o-visualization ancombc-health-state-barplot-relaxed.qzv

    # Export raw stats table so q-values/lfc can be inspected directly
    qiime tools export \
        --input-path ancombc-health-state.qza \
        --output-path exported-asv

    # Collapse to genus level (level 6) for a more interpretable comparison
    qiime taxa collapse \
        --i-table ${table_qza} \
        --i-taxonomy ${taxonomy_qza} \
        --p-level 6 \
        --o-collapsed-table table-l6.qza

    qiime composition ancombc \
        --i-table table-l6.qza \
        --m-metadata-file ${params.metadata} \
        --p-formula 'Health_State' \
        --o-differentials ancombc-genus-health-state.qza \
        --verbose

    qiime composition da-barplot \
        --i-data ancombc-genus-health-state.qza \
        --p-significance-threshold 0.05 \
        --o-visualization ancombc-genus-health-state-barplot.qzv

    qiime composition da-barplot \
        --i-data ancombc-genus-health-state.qza \
        --p-significance-threshold 0.2 \
        --o-visualization ancombc-genus-health-state-barplot-relaxed.qzv

    qiime tools export \
        --input-path ancombc-genus-health-state.qza \
        --output-path exported-genus
    """
}
