process CORE_DIVERSITY {
    tag "Core diversity metrics (phylogenetic)"
    publishDir "${params.outdir}/diversity", mode: 'copy'
    input:
    path rooted_tree_qza
    path table_qza
    output:
    path "core-metrics-results/*",                                      emit: all_results
    path "core-metrics-results/faith_pd_vector.qza",                    emit: faith_pd
    path "core-metrics-results/shannon_vector.qza",                     emit: shannon
    path "core-metrics-results/observed_features_vector.qza",           emit: observed_features
    path "core-metrics-results/evenness_vector.qza",                    emit: evenness
    path "core-metrics-results/unweighted_unifrac_distance_matrix.qza", emit: unweighted_unifrac
    path "core-metrics-results/weighted_unifrac_distance_matrix.qza",   emit: weighted_unifrac
    path "core-metrics-results/bray_curtis_distance_matrix.qza",        emit: bray_curtis
    path "core-metrics-results/jaccard_distance_matrix.qza",            emit: jaccard
    path "core-metrics-results/unweighted_unifrac_pcoa_results.qza",    emit: unweighted_unifrac_pcoa
    path "core-metrics-results/weighted_unifrac_pcoa_results.qza",      emit: weighted_unifrac_pcoa
    path "core-metrics-results/unweighted_unifrac_emperor.qzv",         emit: unweighted_unifrac_emperor
    path "core-metrics-results/weighted_unifrac_emperor.qzv",           emit: weighted_unifrac_emperor
    path "core-metrics-results/bray_curtis_emperor.qzv",                emit: bray_curtis_emperor
    path "core-metrics-results/jaccard_emperor.qzv",                    emit: jaccard_emperor
    path "core-metrics-results/rarefied_table.qza",                     emit: rarefied_table
    script:
    """
    export NUMBA_CACHE_DIR=/tmp/numba_cache_\$\$
    export MPLCONFIGDIR=/tmp/mpl_cache_\$\$
    mkdir -p \$NUMBA_CACHE_DIR \$MPLCONFIGDIR
    mkdir -p core-metrics-results

    # 1. Rarefy table
    qiime feature-table rarefy \
        --i-table ${table_qza} \
        --p-sampling-depth ${params.sampling_depth} \
        --o-rarefied-table core-metrics-results/rarefied_table.qza \
        --verbose

    # 2. Alpha diversity metrics
    qiime diversity alpha \
        --i-table core-metrics-results/rarefied_table.qza \
        --p-metric shannon \
        --o-alpha-diversity core-metrics-results/shannon_vector.qza \
        --verbose

    qiime diversity alpha \
        --i-table core-metrics-results/rarefied_table.qza \
        --p-metric observed_features \
        --o-alpha-diversity core-metrics-results/observed_features_vector.qza \
        --verbose

    qiime diversity alpha \
        --i-table core-metrics-results/rarefied_table.qza \
        --p-metric pielou_e \
        --o-alpha-diversity core-metrics-results/evenness_vector.qza \
        --verbose

    qiime diversity alpha-phylogenetic \
        --i-table core-metrics-results/rarefied_table.qza \
        --i-phylogeny ${rooted_tree_qza} \
        --p-metric faith_pd \
        --o-alpha-diversity core-metrics-results/faith_pd_vector.qza \
        --verbose

    # 3. Beta diversity distance matrices
    qiime diversity beta \
        --i-table core-metrics-results/rarefied_table.qza \
        --p-metric braycurtis \
        --o-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
        --verbose

    qiime diversity beta \
        --i-table core-metrics-results/rarefied_table.qza \
        --p-metric jaccard \
        --o-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
        --verbose

    qiime diversity beta-phylogenetic \
        --i-table core-metrics-results/rarefied_table.qza \
        --i-phylogeny ${rooted_tree_qza} \
        --p-metric unweighted_unifrac \
        --p-threads ${task.cpus} \
        --o-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
        --verbose

    qiime diversity beta-phylogenetic \
        --i-table core-metrics-results/rarefied_table.qza \
        --i-phylogeny ${rooted_tree_qza} \
        --p-metric weighted_unifrac \
        --p-threads ${task.cpus} \
        --o-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
        --verbose

    # 4. PCoA
    qiime diversity pcoa \
        --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
        --o-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
        --verbose

    qiime diversity pcoa \
        --i-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
        --o-pcoa core-metrics-results/jaccard_pcoa_results.qza \
        --verbose

    qiime diversity pcoa \
        --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
        --o-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
        --verbose

    qiime diversity pcoa \
        --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
        --o-pcoa core-metrics-results/weighted_unifrac_pcoa_results.qza \
        --verbose

    # 5. Emperor plots
    qiime emperor plot \
        --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
        --m-metadata-file ${params.metadata} \
        --o-visualization core-metrics-results/bray_curtis_emperor.qzv \
        --verbose

    qiime emperor plot \
        --i-pcoa core-metrics-results/jaccard_pcoa_results.qza \
        --m-metadata-file ${params.metadata} \
        --o-visualization core-metrics-results/jaccard_emperor.qzv \
        --verbose

    qiime emperor plot \
        --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
        --m-metadata-file ${params.metadata} \
        --o-visualization core-metrics-results/unweighted_unifrac_emperor.qzv \
        --verbose

    qiime emperor plot \
        --i-pcoa core-metrics-results/weighted_unifrac_pcoa_results.qza \
        --m-metadata-file ${params.metadata} \
        --o-visualization core-metrics-results/weighted_unifrac_emperor.qzv \
        --verbose
    """
}
