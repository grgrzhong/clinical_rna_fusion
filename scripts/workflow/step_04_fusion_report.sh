#!/bin/bash

#############################################################################
# Clinical RNA Fusion Analysis Workflow - Fusion Report Step
# This script generates a fusion report using Arriba and STAR-Fusion results.
#############################################################################

# Iterate over all sample directories in STARFusion directory
fusion_report() {

    local sample="$1"

    echo "$(date +"%F") $(date +"%T") - Processing sample = ${sample}"

    # Arriba and STAR fusion data
    arriba_fusion_file="$STAR_FUSION_DIR/$sample/fusions.tsv"
    star_fusion_file="$STAR_FUSION_DIR/$sample/star-fusion.fusion_predictions.tsv"

    # Output directory
    output_dir="$STAR_FUSION_DIR/$sample/fusion_report"
    mkdir -p "$output_dir"

    # Run fusion_report for the sample
    echo "$(date +"%F") $(date +"%T") - Running fusion_report ..."

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
        >& "${output_dir}/fusion_report.log"
    
    # Convert the json report to summary in xlsx
    echo "$(date +"%F") $(date +"%T") - Converting JSON report to XLSX ..."
    singularity exec \
        --bind "${STAR_FUSION_DIR}:${STAR_FUSION_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/r.sif" \
        Rscript "${MODULE_DIR}/convert_fusion_report_json_to_xlsx.R" \
        --input_json "${output_dir}/fusions.json" \
        --output_xlsx "${output_dir}/${sample}_fusion_summary.xlsx"
}

# Export the function for parallel execution
export -f fusion_report

# Process samples in parallel
samples=$(find "${STAR_FUSION_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n")

echo "$samples" |
    parallel \
        --jobs "$PARALLEL_JOBS" \
        fusion_report {}
