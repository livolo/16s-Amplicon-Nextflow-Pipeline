process PHYLOGENY {
    tag "Phylogenetic tree (MAFFT + FastTree)"
    publishDir "${params.outdir}/phylogeny", mode: 'copy'

    input:
    path rep_seqs_qza

    output:
    path "aligned-rep-seqs.qza",  emit: aligned
    path "masked-aligned-rep-seqs.qza", emit: masked
    path "unrooted-tree.qza",     emit: unrooted_tree
    path "rooted-tree.qza",       emit: rooted_tree

    script:
    """
    export NUMBA_CACHE_DIR=/tmp/numba_cache_\$\$
    export MPLCONFIGDIR=/tmp/mpl_cache_\$\$
    mkdir -p \$NUMBA_CACHE_DIR \$MPLCONFIGDIR

    qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences ${rep_seqs_qza} \
        --p-n-threads ${task.cpus} \
        --o-alignment aligned-rep-seqs.qza \
        --o-masked-alignment masked-aligned-rep-seqs.qza \
        --o-tree unrooted-tree.qza \
        --o-rooted-tree rooted-tree.qza \
        --verbose
    """
}
