process DENOISE {
    tag "DADA2 denoising"
    publishDir "${params.outdir}/denoised", mode: 'copy'
    input:
    path demux_qza
    output:
    path "table.qza",            emit: table
    path "table.qzv",            emit: table_qzv
    path "rep-seqs.qza",         emit: rep_seqs
    path "rep-seqs.qzv",         emit: rep_seqs_qzv
    path "denoising-stats.qza",  emit: stats
    path "denoising-stats.qzv",  emit: stats_qzv
    script:
    """
    export NUMBA_CACHE_DIR=/tmp/numba_cache_\$\$
    export MPLCONFIGDIR=/tmp/mpl_cache_\$\$
    mkdir -p \$NUMBA_CACHE_DIR \$MPLCONFIGDIR

    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs ${demux_qza} \
        --p-trunc-len-f ${params.trunc_len_f} \
        --p-trunc-len-r ${params.trunc_len_r} \
        --p-trim-left-f ${params.trim_left_f} \
        --p-trim-left-r ${params.trim_left_r} \
        --p-n-threads ${task.cpus} \
        --o-table table.qza \
        --o-representative-sequences rep-seqs.qza \
        --o-denoising-stats denoising-stats.qza \
        --verbose

    qiime feature-table summarize \
        --i-table table.qza \
        --o-visualization table.qzv

    qiime feature-table tabulate-seqs \
        --i-data rep-seqs.qza \
        --o-visualization rep-seqs.qzv

    qiime metadata tabulate \
        --m-input-file denoising-stats.qza \
        --o-visualization denoising-stats.qzv
    """
}
