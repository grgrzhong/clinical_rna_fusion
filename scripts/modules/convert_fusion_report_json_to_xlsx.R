#!/usr/bin/env Rscript

###############################################################################
# Collect cohort level fusion report results
# Author: Zhong Guorui, Eiko
# Date: 2025-07-09
# Usage: input_json = Directory containing fusion report results from Arriba and STAR-Fusion
#        output_xlsx = Path to output xlsx file that contains all fusion results
################################################################################

# Load required libraries
suppressPackageStartupMessages(
    suppressWarnings({
        library(optparse)
        library(jsonlite)
        library(fs)
        library(readxl)
        library(writexl)
        library(tidyverse)
    })
)

# Define command line options
option_list <- list(
    make_option(c("--input_json"),
        type = "character",
        default = NULL,
        help = "JSON file contains fusion report results from Arriba and STAR-Fusion",
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
if (is.null(opt$input_json) || is.null(opt$output_xlsx)) {
    print_help(opt_parser)
    stop("Both --input_json and --output_xlsx are required", call. = FALSE)
}

# Check if input json file exists
if (!file.exists(opt$input_json)) {
    stop(paste("Input JSON file not found:", opt$input_json), call. = FALSE)
}

# Load data
data_flatten <- fromJSON(opt$input_json, flatten = TRUE)

# If data_flattent is empty, then save an empty xlsx file
if (length(data_flatten) == 0) {
    
    write_xlsx(data.frame(), opt$output_xlsx)
    
    message("No data found in the input JSON file, creating an empty xlsx file ...")

} else {

    # Fill in missing position fields
    if (!"starfusion.position" %in% names(data_flatten)) data_flatten$starfusion.position <- NA
    if (!"arriba.position" %in% names(data_flatten)) data_flatten$arriba.position <- NA

    # Clean position strings for comparison
    fusion_summary <- data_flatten |>
        rowwise() |>
        mutate(
            cleaned_starfusion = gsub(":-|:\\+", "", starfusion.position),
            cleaned_starfusion = gsub("-#", "#", cleaned_starfusion),
            position_match = (!is.na(starfusion.position) & !is.na(arriba.position)) &&
                cleaned_starfusion == arriba.position,
            Number_of_callers = case_when(
                position_match ~ 2,
                !is.na(starfusion.position) & is.na(arriba.position) ~ 1,
                is.na(starfusion.position) & !is.na(arriba.position) ~ 1,
                TRUE ~ 0
            ),
            Caller_names = paste(
                c(
                    if (!is.na(starfusion.position)) "STAR-fusion",
                    if (!is.na(arriba.position)) "Arriba"
                ),
                collapse = ", "
            )
        ) |>
        ungroup() |>
        relocate(Number_of_callers, Caller_names, .after = Fusion)

    # Replace empty Databases
    if ("Databases" %in% names(fusion_summary)) {
        if ("Databases" %in% names(fusion_summary)) {
            fusion_summary$Databases <- sapply(fusion_summary$Databases, function(x) {
                if (is.null(x) || length(x) == 0 || identical(x, character(0))) {
                    "#N/A"
                } else if (is.character(x)) {
                    paste(x, collapse = "; ")
                } else {
                    as.character(x)
                }
            })
        }
    }

    # Drop empty position columns
    if (all(is.na(fusion_summary$starfusion.position))) fusion_summary$starfusion.position <- NULL
    if (all(is.na(fusion_summary$arriba.position))) fusion_summary$arriba.position <- NULL

    # Get sample name from folder
    sample_id <- basename(dirname(dirname(opt$input_json)))

    # Save individual fusion summary
    write_xlsx(fusion_summary, opt$output_xlsx)
}
