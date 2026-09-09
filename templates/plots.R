# ============================================================
# plots.R
# Automated plotting and table generation
# for the 16S QIIME2 final report
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
})

# ============================================================
# GENERAL THEME
# ============================================================

report_theme <- function() {

  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}


# ============================================================
# SAFE CSV READER
# ============================================================

read_csv_safe <- function(file) {

  if (!file.exists(file)) {
    return(NULL)
  }

  tryCatch(
    read_csv(
      file,
      show_col_types = FALSE
    ),
    error = function(e) {
      message(
        "Could not read: ",
        file,
        " | ",
        e$message
      )
      NULL
    }
  )
}


# ============================================================
# ANCOM-BC
# ============================================================

read_ancombc_result <- function(
    lfc_file,
    qval_file,
    pval_file = NULL
) {

  if (!file.exists(lfc_file)) {
    return(NULL)
  }

  lfc <- read_csv_safe(lfc_file)

  if (is.null(lfc) || nrow(lfc) == 0) {
    return(NULL)
  }

  qval <- NULL

  if (!is.null(qval_file) &&
      file.exists(qval_file)) {

    qval <- read_csv_safe(qval_file)
  }

  if (ncol(lfc) >= 1) {
    names(lfc)[1] <- "Feature"
  }

  if (!is.null(qval) &&
      ncol(qval) >= 1) {

    names(qval)[1] <- "Feature"
  }

  if (!is.null(qval) &&
      "Feature" %in% names(qval)) {

    result <- left_join(
      lfc,
      qval,
      by = "Feature",
      suffix = c("_lfc", "_qval")
    )

  } else {

    result <- lfc
  }

  result
}


find_lfc_column <- function(df) {

  candidates <- c(
    "lfc",
    "LFC",
    "lfc_lfc",
    "log_fold_change",
    "log2FoldChange"
  )

  found <- candidates[
    candidates %in% names(df)
  ]

  if (length(found) > 0) {
    return(found[1])
  }

  numeric_cols <- names(df)[
    vapply(
      df,
      is.numeric,
      logical(1)
    )
  ]

  if (length(numeric_cols) > 0) {
    return(numeric_cols[1])
  }

  NULL
}


find_qval_column <- function(df) {

  candidates <- c(
    "q_val",
    "qval",
    "q_value",
    "q",
    "q_val_qval"
  )

  found <- candidates[
    candidates %in% names(df)
  ]

  if (length(found) > 0) {
    return(found[1])
  }

  numeric_cols <- names(df)[
    vapply(
      df,
      is.numeric,
      logical(1)
    )
  ]

  q_candidates <- numeric_cols[
    grepl(
      "q|adj|fdr",
      numeric_cols,
      ignore.case = TRUE
    )
  ]

  if (length(q_candidates) > 0) {
    return(q_candidates[1])
  }

  NULL
}


plot_ancombc_lfc <- function(
    df,
    title = "ANCOM-BC log-fold change"
) {

  if (is.null(df) ||
      nrow(df) == 0) {

    return(
      ggplot() +
        annotate(
          "text",
          x = 1,
          y = 1,
          label = "No ANCOM-BC results available"
        ) +
        theme_void()
    )
  }

  lfc_col <- find_lfc_column(df)

  if (is.null(lfc_col)) {

    return(
      ggplot() +
        annotate(
          "text",
          x = 1,
          y = 1,
          label = "LFC column not found"
        ) +
        theme_void()
    )
  }

  plot_df <- df %>%
    mutate(
      LFC = as.numeric(
        .data[[lfc_col]]
      )
    ) %>%
    filter(!is.na(LFC)) %>%
    arrange(LFC)

  if (nrow(plot_df) > 30) {

    plot_df <- bind_rows(
      slice_head(
        plot_df,
        n = 15
      ),
      slice_tail(
        plot_df,
        n = 15
      )
    )
  }

  ggplot(
    plot_df,
    aes(
      x = reorder(
        Feature,
        LFC
      ),
      y = LFC
    )
  ) +
    geom_col() +
    coord_flip() +
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    labs(
      title = title,
      x = NULL,
      y = "Log-fold change"
    ) +
    report_theme()
}


summarize_ancombc <- function(df) {

  if (is.null(df) ||
      nrow(df) == 0) {

    return(
      tibble(
        Total_features = 0,
        Increased = 0,
        Decreased = 0,
        Significant = 0
      )
    )
  }

  lfc_col <- find_lfc_column(df)
  qval_col <- find_qval_column(df)

  if (is.null(lfc_col)) {

    return(
      tibble(
        Total_features = nrow(df),
        Increased = NA_integer_,
        Decreased = NA_integer_,
        Significant = NA_integer_
      )
    )
  }

  x <- as.numeric(
    df[[lfc_col]]
  )

  if (!is.null(qval_col)) {

    q <- as.numeric(
      df[[qval_col]]
    )

    significant <- !is.na(q) &
      q < 0.05

  } else {

    significant <- rep(
      NA,
      length(x)
    )
  }

  tibble(
    Total_features = length(x),

    Increased = sum(
      x > 0,
      na.rm = TRUE
    ),

    Decreased = sum(
      x < 0,
      na.rm = TRUE
    ),

    Significant = ifelse(
      all(is.na(significant)),
      NA,
      sum(
        significant,
        na.rm = TRUE
      )
    )
  )
}


read_ancombc_directory <- function(
    directory
) {

  if (!dir.exists(directory)) {
    return(NULL)
  }

  lfc <- file.path(
    directory,
    "lfc_slice.csv"
  )

  qval <- file.path(
    directory,
    "q_val_slice.csv"
  )

  pval <- file.path(
    directory,
    "p_val_slice.csv"
  )

  if (!file.exists(lfc)) {
    return(NULL)
  }

  read_ancombc_result(
    lfc_file = lfc,
    qval_file = qval,
    pval_file = pval
  )
}


# ============================================================
# SAVE PLOT
# ============================================================

save_report_plot <- function(
    plot,
    file,
    width = 9,
    height = 6
) {

  ggsave(
    filename = file,
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )

  message(
    "Created plot: ",
    file
  )

  invisible(file)
}


# ============================================================
# QIIME2 PCoA
# ============================================================

extract_qiime2_ordination <- function(qza_file) {

  if (!file.exists(qza_file)) {

    message(
      "PCoA file not found: ",
      qza_file
    )

    return(NULL)
  }

  tmp_dir <- tempfile(
    "qiime2_ordination_"
  )

  dir.create(
    tmp_dir,
    recursive = TRUE
  )

  utils::unzip(
    qza_file,
    files = NULL,
    exdir = tmp_dir
  )

  ordination_files <- list.files(
    tmp_dir,
    pattern = "^ordination\\.txt$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(ordination_files) == 0) {

    message(
      "ordination.txt not found in: ",
      qza_file
    )

    unlink(
      tmp_dir,
      recursive = TRUE
    )

    return(NULL)
  }

  ordination_file <- ordination_files[1]

  lines <- readLines(
    ordination_file,
    warn = FALSE
  )

  # ----------------------------------------------------------
  # Locate Site section
  # ----------------------------------------------------------

  site_start <- which(
    grepl(
      "^Site\\t",
      lines
    )
  )

  if (length(site_start) == 0) {

    message(
      "Site section not found in: ",
      qza_file
    )

    unlink(
      tmp_dir,
      recursive = TRUE
    )

    return(NULL)
  }

  site_start <- site_start[1]

  # ----------------------------------------------------------
  # Locate end of Site section
  # ----------------------------------------------------------

  next_sections <- which(
    seq_along(lines) > site_start &
      grepl(
        "^(Biplot|Site constraints|Species|$)",
        lines
      )
  )

  if (length(next_sections) > 0) {

    site_end <- next_sections[1] - 1

  } else {

    site_end <- length(lines)
  }

  site_lines <- lines[
    site_start:site_end
  ]

  site_lines <- site_lines[
    nzchar(trimws(site_lines))
  ]

  # Remove header
  if (length(site_lines) > 0) {
    site_lines <- site_lines[-1]
  }

  if (length(site_lines) == 0) {

    unlink(
      tmp_dir,
      recursive = TRUE
    )

    return(NULL)
  }

  # ----------------------------------------------------------
  # Read sample coordinates
  # ----------------------------------------------------------

  site_text <- paste(
    site_lines,
    collapse = "\n"
  )

  site_df <- read.delim(
    text = site_text,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (ncol(site_df) < 3) {

    message(
      "Insufficient PCoA coordinate columns in: ",
      qza_file
    )

    unlink(
      tmp_dir,
      recursive = TRUE
    )

    return(NULL)
  }

  # First column = sample ID
  sample_id <- as.character(
    site_df[[1]]
  )

  # Remaining columns = PCoA axes
  coordinate_data <- site_df[
    ,
    -1,
    drop = FALSE
  ]

  coordinate_data <- as.data.frame(
    lapply(
      coordinate_data,
      as.numeric
    )
  )

  names(coordinate_data) <- paste0(
    "PC",
    seq_len(
      ncol(coordinate_data)
    )
  )

  result <- bind_cols(
    tibble(
      SampleID = sample_id
    ),
    as_tibble(
      coordinate_data
    )
  )

  unlink(
    tmp_dir,
    recursive = TRUE
  )

  result
}


# ============================================================
# EXTRACT PCoA VARIANCE
# ============================================================

extract_qiime2_variance <- function(
    qza_file
) {

  if (!file.exists(qza_file)) {
    return(NULL)
  }

  tmp_dir <- tempfile(
    "qiime2_variance_"
  )

  dir.create(
    tmp_dir,
    recursive = TRUE
  )

  utils::unzip(
    qza_file,
    files = NULL,
    exdir = tmp_dir
  )

  ordination_files <- list.files(
    tmp_dir,
    pattern = "^ordination\\.txt$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(ordination_files) == 0) {

    unlink(
      tmp_dir,
      recursive = TRUE
    )

    return(NULL)
  }

  lines <- readLines(
    ordination_files[1],
    warn = FALSE
  )

  prop_start <- which(
    grepl(
      "^Proportion explained\\t",
      lines
    )
  )

  if (length(prop_start) == 0) {

    unlink(
      tmp_dir,
      recursive = TRUE
    )

    return(NULL)
  }

  prop_line <- lines[
    prop_start[1] + 1
  ]

  values <- strsplit(
    prop_line,
    "\t",
    fixed = TRUE
  )[[1]]

  values <- as.numeric(
    values
  )

  unlink(
    tmp_dir,
    recursive = TRUE
  )

  values
}


# ============================================================
# CREATE PCoA PLOT
# ============================================================

plot_qiime2_pcoa <- function(
    qza_file,
    metadata,
    group_column,
    title
) {

  coords <- extract_qiime2_ordination(
    qza_file
  )

  variance <- extract_qiime2_variance(
    qza_file
  )

  if (is.null(coords)) {

    return(
      ggplot() +
        annotate(
          "text",
          x = 1,
          y = 1,
          label = "PCoA coordinates unavailable"
        ) +
        theme_void()
    )
  }

  # ----------------------------------------------------------
  # Prepare metadata
  # ----------------------------------------------------------

  if (
    !is.null(metadata) &&
    group_column %in% names(metadata)
  ) {

    # IMPORTANT:
    # QIIME2 uses "sample-id"
    # as the metadata sample identifier.
    #
    # Convert it to the internal name SampleID
    # so it matches the PCoA coordinates.

    if (!"sample-id" %in% names(metadata)) {

      message(
        "Metadata column 'sample-id' not found."
      )

      coords[[group_column]] <- "All samples"

    } else {

      metadata_join <- metadata %>%
        mutate(
          SampleID = as.character(
            .data[["sample-id"]]
          )
        ) %>%
        select(
          SampleID,
          all_of(group_column)
        )

      coords <- coords %>%
        mutate(
          SampleID = as.character(
            SampleID
          )
        ) %>%
        left_join(
          metadata_join,
          by = "SampleID"
        )

      # Convert missing groups to Unknown
      coords[[group_column]] <- ifelse(
        is.na(coords[[group_column]]) |
          coords[[group_column]] == "",
        "Unknown",
        as.character(
          coords[[group_column]]
        )
      )
    }

  } else {

    coords[[group_column]] <- "All samples"
  }


  # ----------------------------------------------------------
  # Check PC1 and PC2
  # ----------------------------------------------------------

  if (
    !"PC1" %in% names(coords) ||
    !"PC2" %in% names(coords)
  ) {

    return(
      ggplot() +
        annotate(
          "text",
          x = 1,
          y = 1,
          label = "PC1/PC2 coordinates unavailable"
        ) +
        theme_void()
    )
  }


  # ----------------------------------------------------------
  # Variance explained
  # ----------------------------------------------------------

  pc1_percent <- if (
    !is.null(variance) &&
    length(variance) >= 1
  ) {

    round(
      variance[1] * 100,
      2
    )

  } else {

    NA
  }


  pc2_percent <- if (
    !is.null(variance) &&
    length(variance) >= 2
  ) {

    round(
      variance[2] * 100,
      2
    )

  } else {

    NA
  }


  x_label <- if (
    is.na(pc1_percent)
  ) {

    "PC1"

  } else {

    paste0(
      "PC1 (",
      pc1_percent,
      "%)"
    )
  }


  y_label <- if (
    is.na(pc2_percent)
  ) {

    "PC2"

  } else {

    paste0(
      "PC2 (",
      pc2_percent,
      "%)"
    )
  }


  # ----------------------------------------------------------
  # PCoA plot
  # ----------------------------------------------------------

  ggplot(
    coords,
    aes(
      x = PC1,
      y = PC2,
      color = .data[[group_column]]
    )
  ) +

    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +

    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +

    geom_point(
      size = 4,
      alpha = 0.85
    ) +

    geom_text(
      aes(
        label = SampleID
      ),
      vjust = -0.8,
      size = 3,
      show.legend = FALSE
    ) +

    labs(
      title = title,
      x = x_label,
      y = y_label,
      color = group_column
    ) +

    report_theme()
}


# ============================================================
# ANCOM-BC SUMMARY
# ============================================================

make_ancombc_summary <- function(
    asv_results,
    genus_results
) {

  asv_summary <- summarize_ancombc(
    asv_results
  )

  genus_summary <- summarize_ancombc(
    genus_results
  )

  bind_rows(
    ASV = asv_summary,
    Genus = genus_summary,
    .id = "Level"
  )
}


# ============================================================
# COMMAND LINE ARGUMENTS
# ============================================================

args <- commandArgs(
  trailingOnly = TRUE
)


if (length(args) < 4) {

  stop(
    paste(
      "Usage:",
      "Rscript plots.R",
      "<results_dir>",
      "<output_dir>",
      "<metadata>",
      "<group_column>"
    )
  )
}


results_dir <- args[1]
output_dir <- args[2]
metadata_file <- args[3]
group_column <- args[4]


cat("\n")
cat("==================================================\n")
cat("16S REPORT PLOT GENERATION\n")
cat("==================================================\n")

cat(
  "Results directory :",
  results_dir,
  "\n"
)

cat(
  "Output directory  :",
  output_dir,
  "\n"
)

cat(
  "Metadata          :",
  metadata_file,
  "\n"
)

cat(
  "Group column      :",
  group_column,
  "\n"
)


# ============================================================
# OUTPUT DIRECTORIES
# ============================================================

figures_dir <- file.path(
  output_dir,
  "figures"
)

tables_dir <- file.path(
  output_dir,
  "tables"
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# READ METADATA
# ============================================================

metadata <- tryCatch(

  read_tsv(
    metadata_file,
    show_col_types = FALSE
  ),

  error = function(e) {

    message(
      "Could not read metadata: ",
      metadata_file,
      " | ",
      e$message
    )

    NULL
  }
)


if (!is.null(metadata)) {

  cat(
    "Metadata rows:",
    nrow(metadata),
    "\n"
  )

  cat(
    "Metadata columns:",
    ncol(metadata),
    "\n"
  )

  cat(
    "Metadata columns:",
    paste(
      names(metadata),
      collapse = ", "
    ),
    "\n"
  )

  if (group_column %in% names(metadata)) {

    cat(
      "Group column found:",
      group_column,
      "\n"
    )

    group_table <- metadata %>%
      count(
        .data[[group_column]],
        name = "Samples"
      )

    write_csv(
      group_table,
      file.path(
        tables_dir,
        "group_summary.csv"
      )
    )

  } else {

    warning(
      "Group column not found in metadata: ",
      group_column
    )
  }

} else {

  warning(
    "Metadata could not be read."
  )
}


# ============================================================
# ANCOM-BC INPUTS
# ============================================================

asv_dir <- file.path(
  results_dir,
  "differential-abundance",
  "exported-asv"
)

genus_dir <- file.path(
  results_dir,
  "differential-abundance",
  "exported-genus"
)


asv_results <- read_ancombc_directory(
  asv_dir
)

genus_results <- read_ancombc_directory(
  genus_dir
)


cat("\n")
cat("ANCOM-BC ASV results: ")

if (is.null(asv_results)) {

  cat("NOT FOUND\n")

} else {

  cat(
    nrow(asv_results),
    "features\n"
  )
}


cat("ANCOM-BC Genus results: ")

if (is.null(genus_results)) {

  cat("NOT FOUND\n")

} else {

  cat(
    nrow(genus_results),
    "features\n"
  )
}


# ============================================================
# ASV LFC PLOT
# ============================================================

if (!is.null(asv_results)) {

  p_asv <- plot_ancombc_lfc(
    asv_results,
    title =
      "ANCOM-BC Differential Abundance - ASV Level"
  )

  save_report_plot(
    p_asv,
    file.path(
      figures_dir,
      "ancombc_asv_lfc.png"
    )
  )
}


# ============================================================
# GENUS LFC PLOT
# ============================================================

if (!is.null(genus_results)) {

  p_genus <- plot_ancombc_lfc(
    genus_results,
    title =
      "ANCOM-BC Differential Abundance - Genus Level"
  )

  save_report_plot(
    p_genus,
    file.path(
      figures_dir,
      "ancombc_genus_lfc.png"
    )
  )
}


# ============================================================
# ANCOM-BC SUMMARY TABLE
# ============================================================

ancombc_summary <- make_ancombc_summary(
  asv_results,
  genus_results
)

write_csv(
  ancombc_summary,
  file.path(
    tables_dir,
    "ancombc_summary.csv"
  )
)

cat(
  "Created ANCOM-BC summary table\n"
)


# ============================================================
# SAVE COMPLETE ANCOM-BC TABLES
# ============================================================

if (!is.null(asv_results)) {

  write_csv(
    asv_results,
    file.path(
      tables_dir,
      "ancombc_asv_results.csv"
    )
  )
}


if (!is.null(genus_results)) {

  write_csv(
    genus_results,
    file.path(
      tables_dir,
      "ancombc_genus_results.csv"
    )
  )
}


# ============================================================
# PCoA FIGURES
# ============================================================

pcoa_dir <- file.path(
  results_dir,
  "diversity",
  "core-metrics-results"
)


pcoa_files <- list(

  "Bray-Curtis" =
    file.path(
      pcoa_dir,
      "bray_curtis_pcoa_results.qza"
    ),

  "Jaccard" =
    file.path(
      pcoa_dir,
      "jaccard_pcoa_results.qza"
    ),

  "Weighted UniFrac" =
    file.path(
      pcoa_dir,
      "weighted_unifrac_pcoa_results.qza"
    ),

  "Unweighted UniFrac" =
    file.path(
      pcoa_dir,
      "unweighted_unifrac_pcoa_results.qza"
    )
)


pcoa_output_names <- c(

  "Bray-Curtis" =
    "bray_curtis_pcoa.png",

  "Jaccard" =
    "jaccard_pcoa.png",

  "Weighted UniFrac" =
    "weighted_unifrac_pcoa.png",

  "Unweighted UniFrac" =
    "unweighted_unifrac_pcoa.png"
)


for (method in names(pcoa_files)) {

  qza_file <- pcoa_files[[method]]

  if (!file.exists(qza_file)) {

    message(
      "PCoA QZA not found: ",
      qza_file
    )

    next
  }

  p <- plot_qiime2_pcoa(

    qza_file = qza_file,

    metadata = metadata,

    group_column = group_column,

    title = paste(
      method,
      "PCoA"
    )
  )

  save_report_plot(

    p,

    file.path(
      figures_dir,
      pcoa_output_names[[method]]
    ),

    width = 9,

    height = 7
  )
}


# ============================================================
# RESULTS INVENTORY
# ============================================================

result_files <- c(

  "denoised/denoising-stats.qza",

  "denoised/rep-seqs.qza",

  "denoised/table.qza",

  "imported_data/demux.qza",

  "phylogeny/rooted-tree.qza",

  "taxonomy/taxonomy.qza",

  "taxonomy/taxa-bar-plots.qzv",

  "diversity/core-metrics-results/bray_curtis_pcoa_results.qza",

  "diversity/core-metrics-results/jaccard_pcoa_results.qza",

  "diversity/core-metrics-results/weighted_unifrac_pcoa_results.qza",

  "diversity/core-metrics-results/unweighted_unifrac_pcoa_results.qza"
)


inventory <- tibble(

  File = result_files,

  Exists =
    file.exists(
      file.path(
        results_dir,
        result_files
      )
    )
)


write_csv(

  inventory,

  file.path(
    tables_dir,
    "results_inventory.csv"
  )
)


# ============================================================
# FINAL MESSAGE
# ============================================================

cat("\n")
cat("==================================================\n")
cat("PLOT GENERATION COMPLETED\n")
cat("==================================================\n")


cat("\nFigures created:\n")


figure_files <- list.files(
  figures_dir,
  full.names = TRUE
)


if (length(figure_files) == 0) {

  cat(
    "  No figures generated\n"
  )

} else {

  for (f in figure_files) {

    cat(
      "  ",
      f,
      "\n",
      sep = ""
    )
  }
}


cat("\nTables created:\n")


table_files <- list.files(
  tables_dir,
  full.names = TRUE
)


if (length(table_files) == 0) {

  cat(
    "  No tables generated\n"
  )

} else {

  for (f in table_files) {

    cat(
      "  ",
      f,
      "\n",
      sep = ""
    )
  }
}


cat("\n")
cat("==================================================\n")
