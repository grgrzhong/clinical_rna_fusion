#!/bin/bash
#SBATCH --job-name=clinical_rna_Fusion
#SBATCH --partition=amd
#SBATCH --time=48:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/home/zhonggr/projects/clinical_rna_fusion/slurm/%x_%j.out
#SBATCH --error=/home/zhonggr/projects/clinical_rna_fusion/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

# Load configuration
source "$(dirname "${BASH_SOURCE[0]}")/conf/config.sh"

# Add better error handling
set -e  # Exit on any error
set -u  # Exit on undefined variables
set -o pipefail  # Exit if any command in pipeline fails

echo "======================================================================="
echo "Clinical RNA Fusion Analysis Workflow - Pipeline Execution"
echo "======================================================================="

# Log the time of start
start_time=$(date +"%F %T")
echo "Pipeline started at: $start_time"

# Step 1: Preprocessing
echo "Step 1: Preprocessing and QC..."
bash "${PROJECT_DIR}/scripts/workflow/step_01_preprocess.sh"
if [ $? -ne 0 ]; then
    echo "✗ Error: Preprocessing steps failed. Exiting pipeline."
    exit 1
fi
echo "✓ Preprocessing steps completed"

# Step 2: Arriba fusion detection
echo "Step 2: Arriba fusion detection..."
bash "${PROJECT_DIR}/scripts/workflow/step_02_arriba_fusion.sh"
if [ $? -ne 0 ]; then
    echo "✗ Error: Arriba fusion steps failed. Exiting pipeline."
    exit 1
fi
echo "✓ Arriba fusion steps completed"

# Step 3: STAR-Fusion detection
echo "Step 3: STAR-Fusion detection..."
bash "${PROJECT_DIR}/scripts/workflow/step_03_star_fusion.sh"
if [ $? -ne 0 ]; then
    echo "✗ Error: STAR fusion steps failed. Exiting pipeline."
    exit 1
fi
echo "✓ STAR fusion steps completed"

# Step 4: Generate Fusion report
echo "Step 4: Generating fusion report..."
bash "${PROJECT_DIR}/scripts/workflow/step_04_fusion_report.sh"
if [ $? -ne 0 ]; then
    echo "✗ Error: Fusion report generation failed. Exiting pipeline."
    exit 1
fi
echo "✓ Fusion report completed"

# Step 5: QC metrics for HS and RnaSeq
echo "Step 5: Generating QC metrics..."
bash "${PROJECT_DIR}/scripts/workflow/step_05_qc_metrics.sh"
if [ $? -ne 0 ]; then
    echo "✗ Error: QC metrics generation failed. Exiting pipeline."
    exit 1
fi
echo "✓ QC metrics generation completed"

# Step 6: Feature count (optional)
# echo "Step 6: Feature counting..."
# bash "${PROJECT_DIR}/scripts/workflow/step_06_feature_counts.sh"
# echo "✓ Feature counting completed"

# log the time of completion
end_time=$(date +"%F %T")
echo "Pipeline completed at: $end_time"

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
echo "Fusion results directory: $OUTPUT_DIR"