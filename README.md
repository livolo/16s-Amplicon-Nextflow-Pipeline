# 16S Amplicon Nextflow Pipeline

A reproducible **16S rRNA amplicon sequencing analysis pipeline** implemented using **Nextflow** and **QIIME 2**.

The pipeline automates quality control, read preprocessing, ASV generation, phylogenetic analysis, taxonomic classification, diversity analysis, differential abundance analysis, and automated R Markdown reporting.

---

## Overview

This workflow is designed for bacterial 16S rRNA amplicon sequencing data and provides an end-to-end analysis framework from raw paired-end FASTQ files to an automated analysis report.

The workflow was developed using a human gut microbiome 16S dataset from NCBI BioProject **PRJNA601994**. The dataset contains paired-end amplicon sequencing reads generated using an Illumina MiSeq platform. The study used the **V4 variable region of the 16S rRNA gene**.

### Study dataset

| Property            | Description                  |
| ------------------- | ---------------------------- |
| BioProject          | PRJNA601994                  |
| Assay type          | 16S amplicon                 |
| Sequencing platform | Illumina MiSeq               |
| Library layout      | Paired-end                   |
| Sample type         | Stool                        |
| Target              | Bacterial 16S rRNA V4 region |
| Organism            | Human gut metagenome         |
| SRA Study           | SRP242624                    |

The original study included **840 individuals**, consisting of 524 Parkinson's disease samples and 316 healthy controls.

---

## Workflow

```text
Raw paired-end FASTQ
        │
        ▼
      FastQC
        │
        ▼
      fastp
        │
        ▼
   Quality Control
        │
        ▼
   QIIME 2 Import
        │
        ▼
      DADA2
        │
        ├── Feature table
        ├── Representative sequences
        └── Denoising statistics
        │
        ▼
 Phylogenetic Analysis
        │
        ▼
 Taxonomic Classification
        │
        ▼
 Diversity Analysis
        │
        ├── Alpha diversity
        ├── Beta diversity
        ├── Shannon
        ├── Observed Features
        ├── Pielou's Evenness
        ├── Faith's PD
        ├── Bray-Curtis
        ├── Jaccard
        ├── Unweighted UniFrac
        └── Weighted UniFrac
        │
        ▼
 Differential Abundance
        │
        ▼
 Automated R Markdown Report
```

---

## Pipeline Modules

The workflow is organized into independent Nextflow modules.

| Module                           | Function                                                       |
| -------------------------------- | -------------------------------------------------------------- |
| `01_demultiplexing.nf`           | Import/demultiplex sequencing data                             |
| `02_denoising.nf`                | DADA2 denoising and ASV generation                             |
| `03_phylogenetic_tree.nf`        | Multiple sequence alignment and phylogenetic tree construction |
| `04_taxonomy.nf`                 | Taxonomic classification                                       |
| `05_diversity.nf`                | Core diversity analysis                                        |
| `06_alpha_group_significance.nf` | Alpha diversity group significance                             |
| `07_alpha_rarefaction.nf`        | Alpha rarefaction analysis                                     |
| `08_differential_abundance.nf`   | Differential abundance analysis                                |
| `10_r_report.nf`                 | Automated R Markdown report generation                         |
| `fastqc.nf`                      | FASTQ quality assessment                                       |
| `fastp.nf`                       | Read quality filtering and preprocessing                       |

---

## Quality Control

The pipeline performs quality assessment before and after read preprocessing.

### FastQC

FastQC is used to inspect:

* Per-base sequence quality
* Read length
* Per-sequence quality
* Adapter contamination
* Overrepresented sequences
* GC content

### fastp

`fastp` is used for preprocessing and quality filtering of sequencing reads.

The processed reads can subsequently be inspected again with FastQC to confirm improvement in read quality.

---

## DADA2 Denoising

QIIME 2 DADA2 is used to process paired-end reads and generate **amplicon sequence variants (ASVs)**.

The workflow includes:

```text
Input reads
    ↓
Quality filtering
    ↓
Denoising
    ↓
Paired-end merging
    ↓
Chimera removal
    ↓
ASV feature table
    ↓
Representative sequences
```

For the example dataset, the analysis generated approximately **620 unique ASVs**. The representative sequences had an average length of approximately **252 bp**, consistent with a V4 16S amplicon.

Denoising statistics are retained to evaluate filtering, merging, and non-chimeric read retention.

---

## Phylogenetic Analysis

Representative sequences are aligned and used to construct a phylogenetic tree.

The workflow uses the QIIME 2:

```text
qiime phylogeny align-to-tree-mafft-fasttree
```

The resulting phylogenetic tree is required for phylogenetic diversity metrics such as:

* Faith's Phylogenetic Diversity
* Unweighted UniFrac
* Weighted UniFrac

---

## Taxonomic Classification

ASVs are classified using a pretrained QIIME 2 taxonomic classifier.

The workflow can use reference databases such as **SILVA** or other compatible 16S reference databases.

Taxonomic assignments can be summarized at different taxonomic levels:

```text
Kingdom
  ↓
Phylum
  ↓
Class
  ↓
Order
  ↓
Family
  ↓
Genus
  ↓
Species
```

For V4 16S data, the reference classifier should be compatible with the targeted region.

---

## Diversity Analysis

The pipeline performs both alpha and beta diversity analyses.

### Alpha Diversity

The workflow includes:

* Observed Features
* Shannon diversity
* Pielou's Evenness
* Faith's Phylogenetic Diversity

### Beta Diversity

The workflow includes:

* Bray-Curtis
* Jaccard
* Unweighted UniFrac
* Weighted UniFrac

These analyses allow microbial community diversity and composition to be investigated across samples and metadata-defined groups.

---

## Differential Abundance

The pipeline includes a module for downstream differential abundance analysis.

This can be used to investigate taxa whose abundance differs between predefined sample groups when the study design and sample size support statistical comparison.

---

## Automated R Markdown Report

The final stage of the workflow generates an automated report using **R Markdown**.

The report integrates analysis outputs and visualizations into a single HTML report.

The report can include:

* Denoising statistics
* Alpha diversity
* Beta diversity
* PCoA visualizations
* Taxonomic composition
* Genus-level abundance
* Diversity group comparisons
* Differential abundance results

The reporting stage is implemented through:

```text
modules/10_r_report.nf
```

with templates stored under:

```text
templates/
```

---

## Repository Structure

```text
16s-Amplicon-Nextflow-Pipeline/
│
├── main.nf
├── nextflow.config
├── R.def
├── .gitignore
│
├── modules/
│   ├── 01_demultiplexing.nf
│   ├── 02_denoising.nf
│   ├── 03_phylogenetic_tree.nf
│   ├── 04_taxonomy.nf
│   ├── 05_diversity.nf
│   ├── 06_alpha_group_significance.nf
│   ├── 07_alpha_rarefaction.nf
│   ├── 08_differential_abundance.nf
│   ├── 10_r_report.nf
│   ├── fastp.nf
│   └── fastqc.nf
│
└── templates/
    ├── Report.Rmd
    ├── plots.R
    └── test_report.Rmd
```

Large sequencing datasets, intermediate files, generated results, QIIME 2 artifacts, containers, and Nextflow work directories are excluded from version control.

---

## Requirements

The pipeline requires:

* Nextflow
* Singularity/Apptainer
* QIIME 2
* FastQC
* fastp
* R
* R Markdown

Containerized environments can be used to improve reproducibility.

---

## Running the Pipeline

Clone the repository:

```bash
git clone https://github.com/livolo/16s-Amplicon-Nextflow-Pipeline.git
cd 16s-Amplicon-Nextflow-Pipeline
```

Review the configuration:

```bash
cat nextflow.config
```

Run the workflow using the configured Singularity profile:

```bash
nextflow run main.nf -profile singularity
```

For a test/stub execution:

```bash
nextflow run main.nf -profile singularity -stub-run
```

Resume an interrupted or previously completed workflow:

```bash
nextflow run main.nf -profile singularity -resume
```

---

## Input Data

The pipeline is designed for paired-end FASTQ data.

Example naming convention:

```text
sample1_1.fastq.gz
sample1_2.fastq.gz

sample2_1.fastq.gz
sample2_2.fastq.gz
```

Raw sequencing data should be placed in the configured input directory.

Large raw datasets should **not** be committed to GitHub.

---

## Outputs

The workflow generates analysis results including:

```text
results/
├── imported_data/
├── denoised/
├── phylogeny/
├── taxonomy/
├── diversity/
├── differential-abundance/
└── r_report/
```

Typical QIIME 2 outputs include:

* Feature tables
* Representative sequences
* Denoising statistics
* Taxonomy results
* Phylogenetic trees
* Alpha diversity results
* Beta diversity / PCoA results
* Taxonomic visualizations
* Statistical results
* Final HTML report

---

## Reproducibility

The pipeline uses **Nextflow DSL2** to organize individual analysis steps into modular processes.

Containerized software environments are used where applicable, allowing the computational environment to be separated from the host system.

The workflow is designed to support:

```text
Version-controlled pipeline
          +
Containerized software
          +
Configured parameters
          +
Automated execution
          =
Reproducible analysis
```

---

## Important Notes

### Study design

Statistical analyses such as differential abundance, group significance, and PERMANOVA should only be interpreted when the available sample size and experimental design support those comparisons.

### 16S resolution

16S amplicon sequencing provides taxonomic profiling but does not always provide reliable species-level resolution. Species-level assignments should therefore be interpreted according to classifier confidence and the targeted 16S region.

### Generated files

Raw reads, intermediate files, QIIME 2 artifacts, result directories, containers, and Nextflow work directories are intentionally excluded from this repository.

---

## Dataset Reference

The example dataset is associated with:

**NCBI BioProject:** PRJNA601994

**SRA Study:** SRP242624

The dataset metadata identifies the assay as amplicon sequencing, with paired-end Illumina MiSeq reads and stool samples.

---

## Author

**Kirti Vishwakarma**

GitHub: [@livolo](https://github.com/livolo)

---

## License

This repository contains the analysis workflow and code. Dataset ownership and reuse are subject to the terms of the original data source and associated study.
