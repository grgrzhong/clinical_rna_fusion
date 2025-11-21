#!/usr/bin/env Rscript

###############################################################################
## Extract QC metrics from MultiQC output files
## Authors: Zhong Guorui
## Usage: Rscript extract_multiqc_metrics.R <multiqc_data_dir> <output_csv> [sample_name]
################################################################################


# Load required libraries
suppressPackageStartupMessages(
    suppressWarnings({
        library(optparse)
        library(tidyverse)
    })
)

# Define command line options
option_list <- list(
    make_option(c("--hs_metrics"),
        type = "character",
        default = NULL,
        help = "Path to MultiQC Picard HsMetrics file",
        metavar = "character"
    ),
    make_option(c("--rna_metrics"),
        type = "character",
        default = NULL,
        help = "Path to MultiQC Picard RnaSeqMetrics file",
        metavar = "character"
    ),
    make_option(c("--sample_name"),
        type = "character",
        default = NULL,
        help = "Sample name for the output",
        metavar = "character"
    ),
    make_option(c("--output_csv"),
        type = "character",
        default = NULL,
        help = "Output CSV file path",
        metavar = "character"
    )
)

# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$hs_metrics) || is.null(opt$rna_metrics) || is.null(opt$sample_name) || is.null(opt$output_csv)) {
    print_help(opt_parser)
    stop("All arguments are required: --hs_metrics, --rna_metrics, --sample_name, --output_csv", call. = FALSE)
}

# Check if input files exist
if (!file.exists(opt$hs_metrics)) {
    stop(paste("HsMetrics file not found:", opt$hs_metrics), call. = FALSE)
}

if (!file.exists(opt$rna_metrics)) {
    stop(paste("RnaSeqMetrics file not found:", opt$rna_metrics), call. = FALSE)
}

# Define required columns
required_columns <- c(
    "TOTAL_READS",
    "PCT_PF_UQ_READS_ALIGNED",
    "PCT_SELECTED_BASES",
    "PCT_OFF_BAIT",
    "MEAN_TARGET_COVERAGE",
    "PCT_RIBOSOMAL_BASES",
    "PCT_CODING_BASES",
    "PCT_UTR_BASES",
    "PCT_INTRONIC_BASES",
    "PCT_INTERGENIC_BASES",
    "PCT_MRNA_BASES"
)

# Read HsMetrics data
hs_metrics_data <- read.table(
    opt$hs_metrics, 
    header = TRUE, 
    sep = "\t", 
    stringsAsFactors = FALSE, 
    comment.char = ""
) |> as_tibble()

# Extract available HsMetrics columns
hs_metrics_cols <- colnames(hs_metrics_data)
hs_metrics_in <- hs_metrics_cols[hs_metrics_cols %in% required_columns]
hs_metrics_tbl <- hs_metrics_data[, hs_metrics_in, drop = FALSE]

# Read RnaSeqMetrics data
rna_metrics_data <- read.table(
    opt$rna_metrics, 
    header = TRUE, 
    sep = "\t", 
    stringsAsFactors = FALSE, 
    comment.char = ""
) |> as_tibble()

# Extract available RnaSeqMetrics columns
rna_metrics_cols <- colnames(rna_metrics_data)
rna_metrics_in <- rna_metrics_cols[rna_metrics_cols %in% required_columns]
rna_metrics_tbl <- rna_metrics_data[, rna_metrics_in, drop = FALSE]

# Combine HS and RNA metrics
reformated_columns <- c(
    "PCT_PF_UQ_READS_ALIGNED",
    "PCT_SELECTED_BASES",
    "PCT_OFF_BAIT",
    "PCT_RIBOSOMAL_BASES"
)

combined_metrics <- tibble(
    Sample = opt$sample_name
) |>
    bind_cols(
        bind_cols(hs_metrics_tbl, rna_metrics_tbl) |>
            select(all_of(required_columns))
    ) |> 
    mutate(
        across(
            all_of(reformated_columns),
            ~ paste0(round(.x * 100, 3), "%")
        ),
        across(
            matches("PCT") & !all_of(reformated_columns),
            ~ paste0(round(.x, 3), "%")
        )
    )

# Write combined metrics to CSV
write_csv(combined_metrics, opt$output_csv)
