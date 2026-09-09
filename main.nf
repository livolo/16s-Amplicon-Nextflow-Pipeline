/*
 * ============================================================
 * 16S rRNA Amplicon Analysis Pipeline
 * QIIME2 + Nextflow DSL2
 *
 * Workflow:
 *
 * 01  Import paired-end FASTQ
 * 02  DADA2 denoising
 * 03  Phylogenetic tree
 * 04  Taxonomy classification
 * 05  Core diversity metrics
 * 06  Alpha diversity group significance
 * 07  Alpha rarefaction
 * 08  Differential abundance
 * 10  Automated R Markdown final report
 *
 * ============================================================
 */


/*
 * ============================================================
 * MODULE INCLUDES
 * ============================================================
 */

include { IMPORT_DATA } from './modules/01_demultiplexing.nf'

include { DENOISE } from './modules/02_denoising.nf'

include { PHYLOGENY } from './modules/03_phylogenetic_tree.nf'

include { TAXONOMY } from './modules/04_taxonomy.nf'

include { CORE_DIVERSITY } from './modules/05_diversity.nf'

include { ALPHA_GROUP_SIGNIFICANCE } from './modules/06_alpha_group_significance.nf'

include { ALPHA_RAREFACTION } from './modules/07_alpha_rarefaction.nf'

include { DIFFERENTIAL_ABUNDANCE } from './modules/08_differential_abundance.nf'

include { R_REPORT } from './modules/10_r_report.nf'


/*
 * ============================================================
 * WORKFLOW
 * ============================================================
 */

workflow {

    /*
     * --------------------------------------------------------
     * Input channels
     * --------------------------------------------------------
     */

    manifest_ch = Channel.fromPath(
        params.manifest,
        checkIfExists: true
    )

    classifier_ch = Channel.fromPath(
        params.classifier,
        checkIfExists: true
    )

    metadata_ch = Channel.fromPath(
        params.metadata,
        checkIfExists: true
    )


    /*
     * --------------------------------------------------------
     * 01. IMPORT DATA
     * --------------------------------------------------------
     */

    IMPORT_DATA(
        manifest_ch
    )


    /*
     * --------------------------------------------------------
     * 02. DADA2 DENOISING
     * --------------------------------------------------------
     */

    DENOISE(
        IMPORT_DATA.out.qza
    )


    /*
     * --------------------------------------------------------
     * 03. PHYLOGENETIC TREE
     * --------------------------------------------------------
     */

    PHYLOGENY(
        DENOISE.out.rep_seqs
    )


    /*
     * --------------------------------------------------------
     * 04. TAXONOMY CLASSIFICATION
     * --------------------------------------------------------
     */

    TAXONOMY(
        DENOISE.out.rep_seqs,
        classifier_ch,
        DENOISE.out.table
    )


    /*
     * --------------------------------------------------------
     * 05. CORE DIVERSITY
     * --------------------------------------------------------
     */

    CORE_DIVERSITY(
        PHYLOGENY.out.rooted_tree,
        DENOISE.out.table
    )


    /*
     * --------------------------------------------------------
     * 06. ALPHA DIVERSITY GROUP SIGNIFICANCE
     * --------------------------------------------------------
     */

    ALPHA_GROUP_SIGNIFICANCE(
        CORE_DIVERSITY.out.faith_pd,
        CORE_DIVERSITY.out.shannon,
        CORE_DIVERSITY.out.evenness,
        CORE_DIVERSITY.out.observed_features
    )


    /*
     * --------------------------------------------------------
     * 07. ALPHA RAREFACTION
     * --------------------------------------------------------
     */

    ALPHA_RAREFACTION(
        DENOISE.out.table,
        PHYLOGENY.out.rooted_tree
    )


    /*
     * --------------------------------------------------------
     * 08. DIFFERENTIAL ABUNDANCE
     * --------------------------------------------------------
     */

    DIFFERENTIAL_ABUNDANCE(
        DENOISE.out.table,
        TAXONOMY.out.taxonomy
    )


    /*
     * ========================================================
     * 10. FINAL R MARKDOWN REPORT
     *
     * IMPORTANT:
     *
     * The report must NOT start simply because the results
     * directory exists.
     *
     * It waits for outputs from the previous analysis modules.
     * These outputs act as completion signals.
     * ========================================================
     */

    results_ch = Channel.value(
        file("${projectDir}/results")
    )


    R_REPORT(
        results_ch,
        metadata_ch,

        /*
         * Completion signal 1
         * Import
         */
        IMPORT_DATA.out.qza,

        /*
         * Completion signal 2
         * DADA2
         */
        DENOISE.out.table,

        /*
         * Completion signal 3
         * Phylogeny
         */
        PHYLOGENY.out.rooted_tree,

        /*
         * Completion signal 4
         * Taxonomy
         */
        TAXONOMY.out.taxonomy,

        /*
         * Completion signal 5
         * Core diversity
         */
        CORE_DIVERSITY.out.faith_pd,

        /*
         * Completion signal 6
         * Alpha group significance
         */
        ALPHA_GROUP_SIGNIFICANCE.out.faith_pd_sig,

        /*
         * Completion signal 7
         * Alpha rarefaction
         */
        ALPHA_RAREFACTION.out.rarefaction_qzv,

        /*
         * Completion signal 8
         * Differential abundance
         */
        DIFFERENTIAL_ABUNDANCE.out.ancombc_results
    )


    /*
     * --------------------------------------------------------
     * Report output
     * --------------------------------------------------------
     */

    R_REPORT.out.report_html.view {
        "R report: $it"
    }

}


/*
 * ============================================================
 * END OF WORKFLOW
 * ============================================================
 */
