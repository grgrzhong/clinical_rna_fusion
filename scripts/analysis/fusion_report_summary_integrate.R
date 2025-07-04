# Load required libraries
library(jsonlite)
library(dplyr)
library(openxlsx)
library(stringr)

# Define base directory
base_dir <- "/media/hkusarc/Data/STUMP_RNA_seq_ALL/Output_StarFusion"

# Find all fusions.json files
json_files <- list.files(base_dir, pattern = "fusions.json", recursive = TRUE, full.names = TRUE)

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
  fusion_summary <- data_flatten %>%
    rowwise() %>%
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
    ) %>%
    ungroup() %>%
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
combined_summary <- combined_all %>%
  group_by(Fusion, starfusion.position, arriba.position) %>%
  summarise(
    Number_of_callers = max(Number_of_callers),
    Caller_names = paste(unique(unlist(strsplit(paste(Caller_names, collapse = ", "), ", "))), collapse = ", "),
    Detected_in_n_samples = n_distinct(Sample),
    Detected_in_samples = paste(sort(unique(Sample)), collapse = ", "),
    .groups = "drop"
  ) %>%
  relocate(starfusion.position, arriba.position, .after = Caller_names)



# Save master Excel
write.xlsx(combined_summary, file.path(base_dir, "Recurrent_Fusion_Summary.xlsx"), rowNames = FALSE)
cat("Master summary saved to: Recurrent_Fusion_Summary.xlsx\n")
