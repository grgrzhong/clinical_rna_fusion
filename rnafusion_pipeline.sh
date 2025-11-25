#!/bin/bash
#SBATCH --job-name=clinical_rna_fusion
#SBATCH --partition=amd
#SBATCH --time=48:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/lustre1/g/path_my/pipeline/clinical_rna_fusion/slurm/%x_%j.out
#SBATCH --error=/lustre1/g/path_my/pipeline/clinical_rna_fusion/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

# Determine the script directory - works both locally and on SLURM
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    # Running on SLURM - use the submit directory
    PIPELINE_DIR="${SLURM_SUBMIT_DIR}"
else
    # Running locally - use the directory containing this script
    PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

# Load configuration
source "${PIPELINE_DIR}/conf/config.sh" "$@"

# Add better error handling
set -e  # Exit on any error
set -u  # Exit on undefined variables
set -o pipefail  # Exit if any command in pipeline fails

# Log the start time to calculate total time taken
start_time=$(date +"%F %T")

echo "======================================================================="
echo "Clinical RNA Fusion Analysis Workflow - Pipeline"
echo "======================================================================="

# ## Step 1: Preprocessing
# echo "$(date +"%F") $(date +"%T") Step 1: Preprocessing and QC ..."
# bash "${PROJECT_DIR}/scripts/workflow/step_01_preprocess.sh"
# if [ $? -ne 0 ]; then
#     echo "✗ Error: Preprocessing steps failed. Exiting pipeline."
#     exit 1
# fi
# echo "$(date +"%F") $(date +"%T") Step 1: Preprocessing and QC (✓) "

# ## Step 2: Arriba fusion detection
# echo "$(date +"%F") $(date +"%T") Step 2: Arriba fusion detection ..."
# bash "${PROJECT_DIR}/scripts/workflow/step_02_arriba_fusion.sh"
# if [ $? -ne 0 ]; then
#     echo "✗ Error: Arriba fusion steps failed. Exiting pipeline."
#     exit 1
# fi
# echo "$(date +"%F") $(date +"%T") Step 2: Arriba fusion detection (✓)"

# ## Step 3: STAR-Fusion detection
# echo "$(date +"%F") $(date +"%T") Step 3: STAR-Fusion detection ..."
# bash "${PROJECT_DIR}/scripts/workflow/step_03_star_fusion.sh"
# if [ $? -ne 0 ]; then
#     echo "✗ Error: STAR fusion steps failed. Exiting pipeline."
#     exit 1
# fi
# echo "$(date +"%F") $(date +"%T") Step 3: STAR-Fusion detection (✓)"

# ## Step 4: Generate Fusion report
# echo "$(date +"%F") $(date +"%T") Step 4: Generating fusion report ..."
# bash "${PROJECT_DIR}/scripts/workflow/step_04_fusion_report.sh"
# if [ $? -ne 0 ]; then
#     echo "✗ Error: Fusion report generation failed. Exiting pipeline."
#     exit 1
# fi
# echo "$(date +"%F") $(date +"%T") Step 4: Generating fusion report (✓)"

# ## Step 5: QC metrics for HS and RnaSeq
# echo "$(date +"%F") $(date +"%T") Step 5: Generating QC metrics ..."
# bash "${PROJECT_DIR}/scripts/workflow/step_05_qc_metrics.sh"
# if [ $? -ne 0 ]; then
#     echo "✗ Error: QC metrics generation failed. Exiting pipeline."
#     exit 1
# fi
# echo "$(date +"%F") $(date +"%T") Step 5: Generating QC metrics (✓)"

## Step 6: Export QC metrics and fusion reports
echo "$(date +"%F") $(date +"%T") Step 6: Exporting QC metrics and fusion reports ..."
bash "${PROJECT_DIR}/scripts/workflow/step_06_export_reports.sh"
if [ $? -ne 0 ]; then
    echo "✗ Error: QC metrics export failed. Exiting pipeline."
    exit 1
fi
echo "$(date +"%F") $(date +"%T") Step 6: Exporting QC metrics and fusion reports (✓)"

## Step 7: Generate feature counts for expression matrix
echo "$(date +"%F") $(date +"%T") Step 7: Generating feature counts for expression matrix ..."
bash "${PROJECT_DIR}/scripts/workflow/step_07_feature_counts.sh"
if [ $? -ne 0 ]; then
    echo "✗ Error: Feature counts generation failed. Exiting pipeline."
    exit 1
fi
echo "$(date +"%F") $(date +"%T") Step 7: Generating feature counts for expression matrix (✓)"

# log the time of completion
end_time=$(date +"%F %T")
echo "$(date +"%F") $(date +"%T") Pipeline completed at: $end_time"

# Calculate total time taken
start_time_sec=$(date -d "$start_time" +%s)
end_time_sec=$(date -d "$end_time" +%s)
total_time=$((end_time_sec - start_time_sec))
hours=$((total_time / 3600))
minutes=$(((total_time % 3600) / 60))
seconds=$((total_time % 60))

# Print the pipeline summary
echo "======================================================================="
echo "Clinical RNA Fusion Analysis Workflow - Running Summary"
echo "======================================================================="
echo "Start Time:               $start_time"
echo "End Time:                 $end_time"
echo "Total Time:               ${hours}h ${minutes}m ${seconds}s"
echo "Output directory:         $OUTPUT_DIR"