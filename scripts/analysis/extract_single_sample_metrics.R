#!/usr/bin/env Rscript

# R script to extract metrics from a single sample's HS and RNA metrics files
# Usage: Rscript extract_single_sample_metrics.R <sample_name> <hs_metrics_file> <rna_metrics_file> <output_csv>

# Function to extract HS metrics
extract_hs_metrics <- function(file_path) {
  tryCatch({
    if (!file.exists(file_path)) {
      cat("Warning: HS metrics file not found:", file_path, "\n")
      return(NULL)
    }
    
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
  }, error = function(e) {
    cat("Error reading HS metrics file", file_path, ":", e$message, "\n")
    return(NULL)
  })
}

# Function to extract RNA metrics
extract_rna_metrics <- function(file_path) {
  tryCatch({
    if (!file.exists(file_path)) {
      cat("Warning: RNA metrics file not found:", file_path, "\n")
      return(NULL)
    }
    
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
  }, error = function(e) {
    cat("Error reading RNA metrics file", file_path, ":", e$message, "\n")
    return(NULL)
  })
}

# Main function
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) != 4) {
    cat("Usage: Rscript extract_single_sample_metrics.R <sample_name> <hs_metrics_file> <rna_metrics_file> <output_csv>\n")
    cat("  sample_name: Name of the sample\n")
    cat("  hs_metrics_file: Path to the HS metrics file\n")
    cat("  rna_metrics_file: Path to the RNA metrics file\n")
    cat("  output_csv: Path to the output CSV file\n")
    quit(status = 1)
  }
  
  sample_name <- args[1]
  hs_file <- args[2]
  rna_file <- args[3]
  output_csv <- args[4]
  
  cat("Processing sample:", sample_name, "\n")
  cat("HS metrics file:", hs_file, "\n")
  cat("RNA metrics file:", rna_file, "\n")
  cat("Output CSV:", output_csv, "\n")
  
  # Initialize row data
  row_data <- data.frame(Sample = sample_name, stringsAsFactors = FALSE)
  
  # Get HS metrics
  hs_data <- extract_hs_metrics(hs_file)
  if (!is.null(hs_data)) {
    row_data <- cbind(row_data, hs_data)
    cat("✓ HS metrics extracted successfully\n")
  } else {
    cat("✗ Could not extract HS metrics\n")
    # Add empty columns for HS metrics
    row_data$TOTAL_READS <- ""
    row_data$PCT_PF_UQ_READS_ALIGNED <- ""
    row_data$PCT_SELECTED_BASES <- ""
    row_data$PCT_OFF_BAIT <- ""
    row_data$MEAN_TARGET_COVERAGE <- ""
  }
  
  # Get RNA metrics
  rna_data <- extract_rna_metrics(rna_file)
  if (!is.null(rna_data)) {
    row_data <- cbind(row_data, rna_data)
    cat("✓ RNA metrics extracted successfully\n")
  } else {
    cat("✗ Could not extract RNA metrics\n")
    # Add empty columns for RNA metrics
    row_data$PCT_RIBOSOMAL_BASES <- ""
    row_data$PCT_CODING_BASES <- ""
    row_data$PCT_UTR_BASES <- ""
    row_data$PCT_INTRONIC_BASES <- ""
    row_data$PCT_INTERGENIC_BASES <- ""
    row_data$PCT_MRNA_BASES <- ""
  }
  
  # Ensure all required columns are present in the correct order
  required_columns <- c(
    "Sample", "TOTAL_READS", "PCT_PF_UQ_READS_ALIGNED", "PCT_SELECTED_BASES",
    "PCT_OFF_BAIT", "MEAN_TARGET_COVERAGE", "PCT_RIBOSOMAL_BASES",
    "PCT_CODING_BASES", "PCT_UTR_BASES", "PCT_INTRONIC_BASES",
    "PCT_INTERGENIC_BASES", "PCT_MRNA_BASES"
  )
  
  # Add missing columns with empty values
  for (col in required_columns) {
    if (!(col %in% names(row_data))) {
      row_data[[col]] <- ""
    }
  }
  
  # Reorder columns
  row_data <- row_data[, required_columns]
  
  # Create output directory if it doesn't exist
  output_dir <- dirname(output_csv)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Check if output file exists to determine if we need to write header
  write_header <- !file.exists(output_csv)
  
  # Write to CSV (append mode)
  write.table(row_data, 
              file = output_csv, 
              sep = ",", 
              row.names = FALSE, 
              col.names = write_header,
              append = !write_header,
              quote = FALSE)
  
  cat("✓ Metrics saved to:", output_csv, "\n")
  
  # Print summary
  cat("\nSample metrics summary:\n")
  for (col in names(row_data)) {
    if (col != "Sample") {
      cat("  ", col, ":", row_data[[col]], "\n")
    }
  }
}

# Run main function
main()
