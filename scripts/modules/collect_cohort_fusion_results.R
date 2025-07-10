#!/usr/bin/env Rscript

###############################################################################
# Collect cohort level fusion report results
# Author: Zhong Guorui, Eiko
# Date: 2025-07-09
# Usage: input_dir = Directory containing fusion report results from Arriba and STAR-Fusion
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
    make_option(c("--input_dir"),
        type = "character",
        default = NULL,
        help = "Directory contains fusion report results from Arriba and STAR-Fusion",
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

# Find all fusions.json files
json_files <- list.files(
    opt$input_dir, 
    pattern = "fusions.json", 
    recursive = TRUE, 
    full.names = TRUE
)

# Store all data for later combination
all_fusions_list <- list()

# Loop through each sample file
for (json_file in json_files) {

    # Load data
    data_flatten <- fromJSON(json_file, flatten = TRUE)
    if (!is.data.frame(data_flatten)) {
        cat("Skipping file (not a data frame):", json_file, "\n")
        next
    }

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
    sample_id <- basename(dirname(dirname(json_file)))

    # Save individual fusion summary
    output_path <- file.path(dirname(json_file), paste0(sample_id, "_fusion_summary.xlsx"))
    write.xlsx(fusion_summary, output_path, rowNames = FALSE, keepNA = TRUE)
    cat("Saved:", output_path, "\n")

    # Add sample ID for aggregation
    fusion_summary$Sample <- sample_id
    all_fusions_list[[sample_id]] <- fusion_summary
}

# Combine all data across samples
combined_all <- bind_rows(all_fusions_list)

# Generate master summary
combined_summary <- combined_all |>
    group_by(Fusion, starfusion.position, arriba.position) |>
    summarise(
        Number_of_callers = max(Number_of_callers),
        Caller_names = paste(
            unique(unlist(strsplit(paste(Caller_names, collapse = ", "), ", "))), 
            collapse = ", "
        ),
        Detected_in_n_samples = n_distinct(Sample),
        Detected_in_samples = paste(sort(unique(Sample)), collapse = ", "),
        .groups = "drop"
    ) |>
    relocate(starfusion.position, arriba.position, .after = Caller_names)

# Save summary to xlsx
write_xlsx(
    list(
        "fusion_report" = combined_all,
        "fusion_summary" = combined_summary
    ), 
    opt$output_xlsx
)

# # Get the sample names from the fusion files
# sample_names <- fusion_files |> 
#     path_dir() |> 
#     path_dir() |> 
#     path_file()

# # Load the data from each fusion report file
# fusion_tbl <- tibble(Sample = sample_names, fusion_file = fusion_files) |>
#     mutate(
#         fusion_data = map2(
#             fusion_file, Sample,
#             \(f, s) {
#                 cat("Loading fusion report for sample:", s, "\n")
#                 readxl::read_xlsx(f)
#             }
#         )
#     ) |>
#     unnest(fusion_data) |> 
#     select(-fusion_file)

# # Write xlsx
# write_xlsx(fusion_tbl, opt$output_xlsx)
