#!/usr/bin/env Rscript

# R script to combine HS and RNA metrics into Excel file
# Usage: Rscript combine_metrics.R <reports_dir>

# Load required libraries
suppressPackageStartupMessages({
    if (!require("openxlsx", quietly = TRUE)) {
        install.packages("openxlsx", repos = "https://cran.r-project.org")
        library(openxlsx)
    }
    if (!require("readr", quietly = TRUE)) {
        install.packages("readr", repos = "https://cran.r-project.org")
        library(readr)
    }
})

# Function to extract HS metrics
extract_hs_metrics <- function(file_path) {
    tryCatch(
        {
            lines <- readLines(file_path)

            # Find the data line (not starting with # or ##)
            for (line in lines) {
                if (!startsWith(line, "#") && !startsWith(line, "##") && nchar(trimws(line)) > 0) {
                    if (grepl("\t", line)) {
                        values <- strsplit(line, "\t")[[1]]
                        if (length(values) > 30) {
                            return(data.frame(
                                TOTAL_READS = values[22],
                                PCT_PF_UQ_READS_ALIGNED = values[27],
                                PCT_SELECTED_BASES = values[7],
                                PCT_OFF_BAIT = values[8],
                                MEAN_TARGET_COVERAGE = values[28],
                                stringsAsFactors = FALSE
                            ))
                        }
                    }
                }
            }
            return(NULL)
        },
        error = function(e) {
            cat("Error reading HS metrics file", file_path, ":", e$message, "\n")
            return(NULL)
        }
    )
}

# Function to extract RNA metrics
extract_rna_metrics <- function(file_path) {
    tryCatch(
        {
            lines <- readLines(file_path)

            # Find the data line (not starting with # or ##)
            for (line in lines) {
                if (!startsWith(line, "#") && !startsWith(line, "##") && nchar(trimws(line)) > 0) {
                    if (grepl("\t", line)) {
                        values <- strsplit(line, "\t")[[1]]
                        if (length(values) > 20) {
                            return(data.frame(
                                PCT_RIBOSOMAL_BASES = values[16],
                                PCT_CODING_BASES = values[17],
                                PCT_UTR_BASES = values[18],
                                PCT_INTRONIC_BASES = values[19],
                                PCT_INTERGENIC_BASES = values[20],
                                PCT_MRNA_BASES = values[21],
                                stringsAsFactors = FALSE
                            ))
                        }
                    }
                }
            }
            return(NULL)
        },
        error = function(e) {
            cat("Error reading RNA metrics file", file_path, ":", e$message, "\n")
            return(NULL)
        }
    )
}

# Main function
main <- function() {
    args <- commandArgs(trailingOnly = TRUE)

    if (length(args) != 1) {
        cat("Usage: Rscript combine_metrics.R <reports_dir>\n")
        quit(status = 1)
    }

    reports_dir <- args[1]

    # Initialize data frame
    combined_data <- data.frame()

    # Process each sample
    hs_metrics_dir <- file.path(reports_dir, "hs_metrics")
    rna_metrics_dir <- file.path(reports_dir, "rna_metrics")

    if (!dir.exists(hs_metrics_dir) || !dir.exists(rna_metrics_dir)) {
        cat("Error: Metrics directories not found in", reports_dir, "\n")
        quit(status = 1)
    }

    # Get all sample directories
    sample_dirs <- list.dirs(hs_metrics_dir, full.names = FALSE, recursive = FALSE)

    for (sample_name in sample_dirs) {
        cat("Processing sample:", sample_name, "\n")

        # Initialize row data
        row_data <- data.frame(Sample = sample_name, stringsAsFactors = FALSE)

        # Get HS metrics
        hs_file <- file.path(hs_metrics_dir, sample_name, paste0(sample_name, "_hs_metrics.txt"))
        if (file.exists(hs_file)) {
            hs_data <- extract_hs_metrics(hs_file)
            if (!is.null(hs_data)) {
                row_data <- cbind(row_data, hs_data)
            } else {
                cat("Warning: Could not extract HS metrics for", sample_name, "\n")
            }
        } else {
            cat("Warning: HS metrics file not found for", sample_name, "\n")
        }

        # Get RNA metrics
        rna_file <- file.path(rna_metrics_dir, sample_name, paste0(sample_name, "_rna_metrics.txt"))
        if (file.exists(rna_file)) {
            rna_data <- extract_rna_metrics(rna_file)
            if (!is.null(rna_data)) {
                row_data <- cbind(row_data, rna_data)
            } else {
                cat("Warning: Could not extract RNA metrics for", sample_name, "\n")
            }
        } else {
            cat("Warning: RNA metrics file not found for", sample_name, "\n")
        }

        # Add to combined data
        if (nrow(combined_data) == 0) {
            combined_data <- row_data
        } else {
            combined_data <- rbind.fill(combined_data, row_data)
        }
    }

    # Ensure all required columns are present
    required_columns <- c(
        "Sample", "TOTAL_READS", "PCT_PF_UQ_READS_ALIGNED", "PCT_SELECTED_BASES",
        "PCT_OFF_BAIT", "MEAN_TARGET_COVERAGE", "PCT_RIBOSOMAL_BASES",
        "PCT_CODING_BASES", "PCT_UTR_BASES", "PCT_INTRONIC_BASES",
        "PCT_INTERGENIC_BASES", "PCT_MRNA_BASES"
    )

    # Add missing columns with empty values
    for (col in required_columns) {
        if (!(col %in% names(combined_data))) {
            combined_data[[col]] <- ""
        }
    }

    # Reorder columns
    combined_data <- combined_data[, required_columns]

    # Create output directory
    output_dir <- file.path(reports_dir, "combined_metrics")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # Save to Excel
    output_file <- file.path(output_dir, "QC_Metrics_Summary.xlsx")

    # Create workbook
    wb <- createWorkbook()
    addWorksheet(wb, "QC_Metrics")
    writeData(wb, "QC_Metrics", combined_data)

    # Format headers
    headerStyle <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2")
    addStyle(wb, "QC_Metrics", headerStyle, rows = 1, cols = 1:ncol(combined_data))

    # Auto-size columns
    setColWidths(wb, "QC_Metrics", cols = 1:ncol(combined_data), widths = "auto")

    # Save workbook
    saveWorkbook(wb, output_file, overwrite = TRUE)

    cat("Combined metrics saved to:", output_file, "\n")
    cat("Processed", nrow(combined_data), "samples\n")
}

# Helper function for rbind.fill (in case plyr is not available)
rbind.fill <- function(df1, df2) {
    all_cols <- unique(c(names(df1), names(df2)))

    for (col in all_cols) {
        if (!(col %in% names(df1))) {
            df1[[col]] <- NA
        }
        if (!(col %in% names(df2))) {
            df2[[col]] <- NA
        }
    }

    rbind(df1, df2)
}

# Run main function
main()
