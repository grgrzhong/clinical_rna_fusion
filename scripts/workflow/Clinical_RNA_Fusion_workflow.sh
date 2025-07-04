#!/bin/bash

#############################################################################
# Clinical RNA Fusion Analysis Workflow
#############################################################################
# Simple workflow to run all RNA fusion analysis steps
#############################################################################

# Exit on any error
set -e

# Setup environment variables
export PROJECT_DIR="/mnt/f/projects/clinical_rna_fusion"
export CONTAINER_DIR="${PROJECT_DIR}/containers"
export PRIMARY_SEQ_DIR="${PROJECT_DIR}/data/primary_seq"
export FASTQ_TRIM_DIR="${PROJECT_DIR}/data/fastq_trimmed"
export FASTQC_TRIM_DIR="${PROJECT_DIR}/data/fastqc_trimmed"

export STAR_FUSION_DIR="${PROJECT_DIR}/data/starfusion"
export ARRIBA_DIR="${PROJECT_DIR}/data/arriba"

export REFERENCE_DIR="/mnt/f/projects/Reference"
export STAR_INDEX="${REFERENCE_DIR}/Gencode/STAR_index"
export REFERENCE="${REFERENCE_DIR}/Gencode/gencode.hg38.v36.primary_assembly.fa"
export ANNOTATION="${REFERENCE_DIR}/Gencode/gencode.v36.primary_assembly.annotation.gtf"
export CTAT_RESOURCE_LIB="${REFERENCE_DIR}/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/"

# Create directories if they don't exist
mkdir -p "$CONTAINER_DIR" "$PRIMARY_SEQ_DIR" "$FASTQ_TRIM_DIR" "$FASTQC_TRIM_DIR"

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"



echo "=========================================="
echo "Clinical RNA Fusion Analysis Workflow"
echo "=========================================="
echo "Starting analysis at: $(date)"
echo ""

# Step 1: Setup containers
echo "Step 1: Setting up containers..."
bash "$SCRIPT_DIR/Step_00_Setup_containers.sh"
echo "✓ Container setup completed"
echo ""

# Step 2: Preprocessing
echo "Step 2: Preprocessing and QC..."
bash "$SCRIPT_DIR/Step_01_Preprocess.sh"
echo "✓ Preprocessing completed"
echo ""

# Step 3: STAR-Fusion
echo "Step 3: STAR-Fusion analysis..."
bash "$SCRIPT_DIR/Step_02_Star_fusion.sh"
echo "✓ STAR-Fusion analysis completed"
echo ""

# Step 4: Arriba
echo "Step 4: Arriba fusion detection..."
bash "$SCRIPT_DIR/Step_03_Arriba_fusion.sh"
echo "✓ Arriba analysis completed"
echo ""

# Step 5: Fusion report
echo "Step 5: Generating fusion report..."
bash "$SCRIPT_DIR/Step_04_Fusion_report.sh"
echo "✓ Fusion report completed"
echo ""

# Step 6: QC metrics
echo "Step 6: Collecting QC metrics..."
bash "$SCRIPT_DIR/Step_05_CollectRnaSeqMetrics.sh"
echo "✓ QC metrics collection completed"
echo ""

echo "=========================================="
echo "Clinical RNA Fusion Analysis COMPLETED"
echo "=========================================="
echo "Analysis completed at: $(date)"
echo ""
echo "Check the following directories for results:"
echo "  - outputs/StarFusion/"
echo "  - outputs/Arriba/"
echo "  - results/"