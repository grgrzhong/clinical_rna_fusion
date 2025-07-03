#!/bin/bash

##### prerequisite #####
# conda activate fusion_report_env
# Python version requirement 3.7

##### test #####
# Run the fusion-report - one sample test
# fusion_report run "RNA195" /media/hkusarc/Data/clinical_rnaseq_june11/fusion_report /media/hkusarc/Data/Reference/fusion_report \
#   --arriba /media/hkusarc/Data/clinical_rnaseq_june11/Output-exome-Arriba/RNA195/fusions.tsv \
#   --starfusion /media/hkusarc/Data/clinical_rnaseq_june11/Output-exome-STARFusion/RNA195/star-fusion.fusion_predictions.tsv \
#   --no-cosmic \
#   --allow-multiple-gene-symbols
  

##### fusion_report #####
# Run the fusion-report - all sample 
# Base directories
BASE_DIR="/media/hkusarc/Data/STUMP_RNA_seq_ALL"
ARRIBA_DIR="$BASE_DIR/Output"
STARFUSION_DIR="$BASE_DIR/Output_StarFusion"
REFERENCE_DIR="/media/hkusarc/Data/Reference/fusion_report"

# Iterate over all sample directories in STARFusion directory
for SAMPLE_DIR in "$STARFUSION_DIR"/*; do
  if [ -d "$SAMPLE_DIR" ]; then
    # Extract sample name dynamically
    SAMPLE_NAME=$(basename "$SAMPLE_DIR")
    
    # Define input files
    ARRIBA_FILE="$ARRIBA_DIR/$SAMPLE_NAME/fusions.tsv"
    STARFUSION_FILE="$STARFUSION_DIR/$SAMPLE_NAME/star-fusion.fusion_predictions.tsv"
    
    # Extract sample names from file paths
    ARRIBA_SAMPLE_NAME=$(basename "$(dirname "$ARRIBA_FILE")")
    STARFUSION_SAMPLE_NAME=$(basename "$(dirname "$STARFUSION_FILE")")
    
    # Validate that the sample name matches between ARRIBA and STARFusion files
    if [[ "$ARRIBA_SAMPLE_NAME" == "$SAMPLE_NAME" && "$STARFUSION_SAMPLE_NAME" == "$SAMPLE_NAME" ]]; then
      # Define output directory for the sample
      SAMPLE_OUTPUT_DIR="$STARFUSION_DIR/$SAMPLE_NAME/fusion_report"
      mkdir -p "$SAMPLE_OUTPUT_DIR"
      
      # Run fusion_report for the sample
      fusion_report run "$SAMPLE_NAME" "$SAMPLE_OUTPUT_DIR" "$REFERENCE_DIR" \
        --arriba "$ARRIBA_FILE" \
        --starfusion "$STARFUSION_FILE" \
        --no-cosmic \
        --allow-multiple-gene-symbols
    else
      echo "Warning: Missing input files for sample $SAMPLE_NAME. Skipping..."
    fi
  fi
done