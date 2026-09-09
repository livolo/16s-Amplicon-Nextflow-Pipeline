process IMPORT_DATA {
    
    publishDir "${params.outdir}/imported_data", mode: 'copy'

    input:
    path manifest

    output:
    path "demux.qza", emit: qza
    path "demux.qzv", emit: qzv

    script:
    """
    export NUMBA_CACHE_DIR=/tmp/numba_cache_\$\$
    export MPLCONFIGDIR=/tmp/mpl_cache_\$\$
    mkdir -p \$NUMBA_CACHE_DIR \$MPLCONFIGDIR

    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path ${manifest} \
        --output-path demux.qza \
        --input-format PairedEndFastqManifestPhred33V2

    qiime demux summarize \
        --i-data demux.qza \
        --o-visualization demux.qzv
    """
}
