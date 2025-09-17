#! /bin/bash

#############################################################################
# Clinical RNA Fusion Analysis Workflow - QC Metrics Step
# This script collects gene fusion and qc metrics for entire cohorts
############################################################################

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJECT_DIR

export CONTAINER_DIR="${PROJECT_DIR}/containers"
export MODULE_DIR="${PROJECT_DIR}/scripts/modules"
export DATA_DIR="${PROJECT_DIR}/results/DFSP"
export REPORTS_DIR="${DATA_DIR}/Reports"

mkdir -p "$REPORTS_DIR"

## Export qc reports for entire cohort
singularity exec \
    --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
    --bind "${DATA_DIR}:${DATA_DIR}" \
    --bind /tmp:/tmp \
    "${CONTAINER_DIR}/r.sif" \
    Rscript "${MODULE_DIR}/collect_cohort_qc_metrics.R" \
        --input_dir "${DATA_DIR}" \
        --output_xlsx "${REPORTS_DIR}/qc_metrics_cohort.xlsx"

## Export fusion results for entire cohort
singularity exec \
    --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
    --bind "${DATA_DIR}:${DATA_DIR}" \
    --bind /tmp:/tmp \
    "${CONTAINER_DIR}/r.sif" \
    Rscript "${MODULE_DIR}/collect_cohort_fusion_reports.R" \
        --input_dir "${DATA_DIR}" \
        --output_xlsx "${REPORTS_DIR}/fusion_report_cohort.xlsx"