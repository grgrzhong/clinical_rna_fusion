# Load required libraries
library(jsonlite)
library(dplyr)
library(openxlsx)

# Define the base directory
base_dir <- "/media/hkusarc/Data/STUMP_RNA_seq_ALL/STAR-Fusion-Singularity-Output"

# Get a list of all fusions.json files in the directory structure
json_files <- list.files(base_dir, pattern = "fusions.json", recursive = TRUE, full.names = TRUE)

# Loop through each file and process it
for (json_file in json_files) {
    # Load the JSON file
    data_flatten <- fromJSON(json_file, flatten = TRUE)

    # Add default columns if starfusion.position or arriba.position do not exist
    if (!"starfusion.position" %in% names(data_flatten)) {
        data_flatten$starfusion.position <- NA
    }
    if (!"arriba.position" %in% names(data_flatten)) {
        data_flatten$arriba.position <- NA
    }

    # Create the fusion_summary data frame
    fusion_summary <- data_flatten %>%
        rowwise() %>%
        mutate(
            # Count the number of callers based on non-NA values in starfusion.position and arriba.position
            Number_of_callers = sum(!is.na(starfusion.position), !is.na(arriba.position)),

            # Display caller names based on the presence of values in starfusion.position and arriba.position
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

    # Replace character(0) with "na" in the Databases column
    if ("Databases" %in% names(fusion_summary)) {
        fusion_summary$Databases <- sapply(fusion_summary$Databases, function(x) {
            if (identical(x, character(0))) {
                "#N/A"
            } else {
                x
            }
        })
    }

    # Remove columns if all rows are NA
    if (all(is.na(fusion_summary$starfusion.position))) {
        fusion_summary$starfusion.position <- NULL
    }
    if (all(is.na(fusion_summary$arriba.position))) {
        fusion_summary$arriba.position <- NULL
    }

    # Extract the subdirectory name (e.g., S13, S14, etc.)
    sub_dir <- basename(dirname(json_file))

    # Define the output file path
    output_file <- file.path(dirname(json_file), paste0(sub_dir, "_fusion_summary.xlsx"))
    print(fusion_summary)
    # Save the result to an Excel file
    write.xlsx(fusion_summary, output_file, rowNames = FALSE, keepNA = TRUE)

    # Log the progress
    cat("Processed and saved:", output_file, "\n")
}

# Optional: Display fusion_caller info only
# fusion_caller <- fusion_summary %>%
#   select(Fusion, Number_of_callers, Caller_names)
# print(fusion_caller)
