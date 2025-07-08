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

#############################################################################
# Clinical RNA Fusion Analysis Workflow
#############################################################################
# Simple workflow to run all RNA fusion analysis steps
#############################################################################
# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate rnafusion

# Exit on any error
set -e

# "=========================================================================="
# Set project directory relative to this script
# "=========================================================================="
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJECT_DIR
export CONTAINER_DIR="${PROJECT_DIR}/containers"
export MODULE_DIR="${PROJECT_DIR}/scripts/modules"

# "=========================================================================="
# Input and output directories
# "=========================================================================="
# Use command line argument for PRIMARY_SEQ_DIR, or default
export PRIMARY_SEQ_DIR="${1:-${PROJECT_DIR}/data/primary_seq}"
export RESULTS_DIR="${2:-${PROJECT_DIR}/results}"
export FASTQ_TRIM_DIR="${RESULTS_DIR}/fastq_trimmed"
export FASTQC_TRIM_DIR="${RESULTS_DIR}/fastqc_trimmed"
export ARRIBA_DIR="${RESULTS_DIR}/output"
export STAR_FUSION_DIR="${RESULTS_DIR}/output"
export FEATURE_COUNTS_DIR="${RESULTS_DIR}/feature_counts"

# Number of threads to use for parallel processing for a sample
export THREADS=8  

# Number of jobs to run in parallel, must be less than the number of samples
export JOBS="${3:-2}"

# Create directories if they don't exist
mkdir -p "$RESULTS_DIR" "$FASTQ_TRIM_DIR" "$FASTQC_TRIM_DIR" 
mkdir -p "$ARRIBA_DIR" "$STAR_FUSION_DIR"

# "=========================================================================="
# Reference and annotation directories
# "=========================================================================="
export REFERENCE_DIR="/mnt/f/projects/Reference"
export STAR_INDEX="${REFERENCE_DIR}/Gencode/STAR_index"
export REFERENCE="${REFERENCE_DIR}/Gencode/gencode.hg38.v36.primary_assembly.fa"
export BAIT_INTERVALS="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2/hg38/xgen-exome-hyb-panel-v2-probes-hg38.interval_list"
export TARGET_INTERVALS="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2/hg38/xgen-exome-hyb-panel-v2-targets-hg38.interval_list"

export REF_FLAT="${REFERENCE_DIR}/Picard_QC/CollectRnaSeqMetrics/refFlat.txt"
export RIBO_INTERVALS="${REFERENCE_DIR}/Picard_QC/CollectRnaSeqMetrics/v44/GRCh38_gencode_v44_rRNA.interval_list"

export ANNOTATION="${REFERENCE_DIR}/Gencode/gencode.v36.primary_assembly.annotation.gtf"
export CTAT_RESOURCE_LIB="${REFERENCE_DIR}/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/"

export BLACKLIST="${REFERENCE_DIR}/arriba_v2.4.0/database/blacklist_hg38_GRCh38_v2.4.0.tsv.gz"
export KNOWN_FUSION="${REFERENCE_DIR}/arriba_v2.4.0/database/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz"
export PROTEIN_DOMAINS="${REFERENCE_DIR}/arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.gff3"
export CYTOBANDS="${REFERENCE_DIR}/arriba_v2.4.0/database/cytobands_hg38_GRCh38_v2.4.0.tsv"

# Print out the environment information
echo "========================================================================"
echo "Clinical RNA Fusion Analysis Workflow"
echo "Project Directory:            $PROJECT_DIR"
echo "Container Directory:          $CONTAINER_DIR"
echo "Module Directory:             $MODULE_DIR"
echo "Primary Seq Directory:        $PRIMARY_SEQ_DIR"
echo "Results Directory:            $RESULTS_DIR"
echo "FastQ Trim Directory:         $FASTQ_TRIM_DIR"
echo "FastQC Trim Directory:        $FASTQC_TRIM_DIR"
echo "Arriba Output Directory:      $ARRIBA_DIR"
echo "STAR-Fusion Output Directory: $STAR_FUSION_DIR"
echo "Feature Counts Directory:     $FEATURE_COUNTS_DIR"
echo "Reference Directory:          $REFERENCE_DIR"
echo "STAR Index:                   $STAR_INDEX"
echo "Reference Genome:             $REFERENCE"
echo "Bait Intervals:               $BAIT_INTERVALS"
echo "Target Intervals:             $TARGET_INTERVALS"
echo "RefFlat:                      $REF_FLAT"
echo "Ribosomal Intervals:          $RIBO_INTERVALS"
echo "Annotation GTF:               $ANNOTATION"
echo "CTAT Resource Library:        $CTAT_RESOURCE_LIB"
echo "Blacklist:                    $BLACKLIST"
echo "Known Fusions:                $KNOWN_FUSION"
echo "Protein Domains:              $PROTEIN_DOMAINS"
echo "Cytobands:                    $CYTOBANDS"
echo "Number of threads:            $THREADS"
echo "Number of GNU parallel jobs:  $JOBS"
echo "======================================================================="
echo "Starting analysis at: $(date)"

# "=========================================================================="
# Main workflow steps
# "=========================================================================="
# Step 0: Setup containers if not already done
if [ ! -d "$CONTAINER_DIR" ]; then
    echo "Container directory not found. Setting up containers..."
    mkdir -p "$CONTAINER_DIR"
    bash "${PROJECT_DIR}/scripts/workflow/step_00_setup_containers.sh"
    echo "✓ Container setup completed"
else
    echo "Container directory already exists: $CONTAINER_DIR"
    echo "Skipping container setup step."
fi

# Step 1: Preprocessing
echo "Step 1: Preprocessing and QC..."
bash "${PROJECT_DIR}/scripts/workflow/step_01_preprocess.sh"
echo "✓ Preprocessing completed"
echo ""

# Step 2: Arriba fusion detection
echo "Step 2: Arriba fusion detection..."
bash "${PROJECT_DIR}/scripts/workflow/step_02_arriba_fusion.sh"
echo "✓ Arriba detection completed"
echo ""

# Step 3: STAR-Fusion detection
echo "Step 3: STAR-Fusion detection..."
bash "${PROJECT_DIR}/scripts/workflow/step_03_star_fusion.sh"
echo "✓ STAR-Fusion detection completed"
echo ""

# Step 4: Generate Fusion report
echo "Step 4: Generating fusion report..."
bash "${PROJECT_DIR}/scripts/workflow/step_04_fusion_report.sh"
echo "✓ Fusion report completed"
echo ""

# Step 5: QC metrics
echo "Step 5: Collecting QC metrics..."
bash "${PROJECT_DIR}/scripts/workflow/step_05_qc_metrics.sh"
echo "✓ QC metrics collection completed"
echo ""

# Step 6: Feature count (optional)
echo "Step 6: Feature counting..."
bash "${PROJECT_DIR}/scripts/workflow/step_06_feature_counts.sh"
echo "✓ Feature counting completed"
echo ""

echo "=========================================="
echo "Clinical RNA Fusion Analysis COMPLETED"
echo "=========================================="
echo "Analysis completed at: $(date)"
echo ""
echo "Check the following directories for results:"
echo "  - ${ARRIBA_DIR}/"
echo "  - ${STAR_FUSION_DIR}/"
echo "  - ${REPORTS_DIR}/"