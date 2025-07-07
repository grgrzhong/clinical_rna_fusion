#!/bin/bash

#############################################################################
# Clinical RNA Fusion Analysis Workflow
#############################################################################
# Simple workflow to run all RNA fusion analysis steps
#############################################################################

# Exit on any error
set -e

# Set script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SCRIPT_DIR

# Input and output directories
export PROJECT_DIR="/mnt/f/projects/clinical_rna_fusion"
export CONTAINER_DIR="${PROJECT_DIR}/containers"
export MODULE_DIR="${PROJECT_DIR}/scripts/modules"
export PRIMARY_SEQ_DIR="${PROJECT_DIR}/data/primary_seq"
export FASTQ_TRIM_DIR="${PROJECT_DIR}/data/fastq_trimmed"
export FASTQC_TRIM_DIR="${PROJECT_DIR}/data/fastqc_trimmed"

export ARRIBA_DIR="${PROJECT_DIR}/data/output"
export STAR_FUSION_DIR="${PROJECT_DIR}/data/output_test"
export REPORTS_DIR="${PROJECT_DIR}/data/reports"

# Reference and annotation directories
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

# Number of threads to use for parallel processing
export THREADS=16  

# Create directories if they don't exist
mkdir -p "$CONTAINER_DIR" "$PRIMARY_SEQ_DIR" "$FASTQ_TRIM_DIR" "$FASTQC_TRIM_DIR" 
mkdir -p "$ARRIBA_DIR" "$STAR_FUSION_DIR"


echo "=========================================="
echo "Clinical RNA Fusion Analysis Workflow"
echo "=========================================="
echo "Starting analysis at: $(date)"
echo ""

# Step 1: Setup containers if not already done
echo "Step 1: Setting up containers..."
bash "$SCRIPT_DIR/step_00_Setup_containers.sh"
echo "✓ Container setup completed"
echo ""

# Step 2: Preprocessing
echo "Step 2: Preprocessing and QC..."
bash "$SCRIPT_DIR/step_01_Preprocess.sh"
echo "✓ Preprocessing completed"
echo ""

# Step 3: Arriba fusion detection
echo "Step 3: Arriba fusion detection ..."
bash "$SCRIPT_DIR/step_03_Arriba_fusion.sh"
echo "✓ Arriba detection completed"
echo ""

# Step 4: STAR-Fusion detection
echo "Step 4: STAR-Fusion detection ..."
bash "$SCRIPT_DIR/step_02_STAR_fusion.sh"
echo "✓ STAR-Fusion detection completed"
echo ""

# Step 5: Generate Fusion report
echo "Step 5: Generating fusion report..."
bash "$SCRIPT_DIR/step_04_Fusion_report.sh"
echo "✓ Fusion report completed"
echo ""

# Step 6: QC metrics
echo "Step 6: Collecting QC metrics..."
bash "$SCRIPT_DIR/step_05_QC.sh"
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