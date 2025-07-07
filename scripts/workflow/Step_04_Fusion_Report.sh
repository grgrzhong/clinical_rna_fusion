#!/bin/bash

#############################################################################
# Clinical RNA Fusion Analysis Workflow - Fusion Report Step
# This script generates a fusion report using Arriba and STAR-Fusion results.
#############################################################################

# Iterate over all sample directories in STARFusion directory
for sample in $(find "$STAR_FUSION_DIR" -mindepth 1 -maxdepth 1 -type d -printf "%f\n"); do

    echo "$(date +"%F") $(date +"%T") - Processing sample = ${sample}"

    # Arriba and STAR fusion data
    arriba_fusion_file="$STAR_FUSION_DIR/$sample/fusions.tsv"
    star_fusion_file="$STAR_FUSION_DIR/$sample/star-fusion.fusion_predictions.tsv"

    # Check if input files exist
    if [[ ! -f "$arriba_fusion_file" ]]; then
        echo "Warning: Arriba fusion file not found: $arriba_fusion_file"
        continue
    fi

    if [[ ! -f "$star_fusion_file" ]]; then
        echo "Warning: STAR-Fusion file not found: $star_fusion_file"
        continue
    fi

    # Output directory
    output_dir="$STAR_FUSION_DIR/$sample/fusion_report"
    mkdir -p "$output_dir"

    # Run fusion_report for the sample
    echo "$(date +"%F") $(date +"%T") - Running fusion_report for sample = ${sample}"

    singularity exec \
        --bind "${STAR_FUSION_DIR}:${STAR_FUSION_DIR}" \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/fusion_report.sif" \
        fusion_report run "$sample" "$output_dir" "$REFERENCE_DIR/fusion_report" \
        --arriba "${arriba_fusion_file}" \
        --starfusion "${star_fusion_file}" \
        --no-cosmic \
        --allow-multiple-gene-symbols \
        >&"${output_dir}/fusion_report.log"
done
