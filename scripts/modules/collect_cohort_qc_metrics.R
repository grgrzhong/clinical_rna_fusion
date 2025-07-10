#!/usr/bin/env Rscript

###############################################################################
# Collect cohort level hs and rna metrics
# Author: Zhong Guorui, Eiko
# Date: 2025-07-09
# Usage: input_dir = Directory containing hs and rna qc metrics
#        output_xlsx = Path to output xlsx file that contains all hs and rna metrics
################################################################################

# Load required libraries
suppressPackageStartupMessages(
    suppressWarnings({
        library(optparse)
        library(fs)
        library(readxl)
        library(writexl)
        library(tidyverse)
    })
)

# Define command line options
option_list <- list(
    make_option(c("--input_dir"),
        type = "character",
        default = NULL,
        help = "Directory contains hs and rna qc metrics",
        metavar = "character"
    ),
    make_option(c("--output_xlsx"),
        type = "character",
        default = NULL,
        help = "Path to output xlsx file",
        metavar = "character"
    )
)

# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$input_dir) || is.null(opt$output_xlsx)) {
    print_help(opt_parser)
    stop("Both --input_dir and --output_xlsx are required", call. = FALSE)
}

# Check if input dir exist
if (!dir.exists(opt$input_dir)) {
    stop(paste("Input directory not found:", opt$input_dir), call. = FALSE)
}

# Find all qc metrics files
metrics_files <- dir_ls(
    opt$input_dir, 
    glob = "*RnaSeqMetrics.csv", 
    recurse = TRUE
)

# Load the data from each fusion report file
metrics_tbl <- tibble(metrics_file = metrics_files) |>
    mutate(
        qc_data = map(
            metrics_file,
            \(f) {
                # cat("Loading qc report for sample:", s, "\n")
                read_csv(f, show_col_types = FALSE)
            }
        )
    ) |>
    unnest(qc_data) |> 
    select(-metrics_file)

# Write xlsx
write_xlsx(metrics_tbl, opt$output_xlsx)
