#!/bin/bash

#############################################################################
# Clinical RNA Fusion Analysis Workflow - Configuration Script
# This script sets up the environment and variables for the RNA fusion analysis
# workflow.
# It configures directories, container paths, and reference files.
#############################################################################

# Activate conda environment
# source $(conda info --base)/etc/profile.d/conda.sh
# conda activate rnafusion

# Exit on any error
set -e

# "=========================================================================="
# Set project directory relative to this script
# "=========================================================================="
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PROJECT_DIR
export CONTAINER_DIR="${PROJECT_DIR}/containers"
export MODULE_DIR="${PROJECT_DIR}/scripts/modules"

# "=========================================================================="
# Input and output directories
# "=========================================================================="
# Use command line argument for PRIMARY_SEQ_DIR, or default
export PRIMARY_SEQ_DIR="${1:-${PROJECT_DIR}/data/Raw}"
export RESULTS_DIR="${2:-${PROJECT_DIR}/data}"
export FASTQ_TRIM_DIR="${RESULTS_DIR}/Input-trimmed"
export FASTQC_TRIM_DIR="${RESULTS_DIR}/FastQC-trimmed"
export ARRIBA_DIR="${RESULTS_DIR}/Output"
export STAR_FUSION_DIR="${RESULTS_DIR}/Output"
export REPORTS_DIR="${RESULTS_DIR}/Reports"
export FEATURE_COUNTS_DIR="${RESULTS_DIR}/Feature-counts"

# Create directories if they don't exist
mkdir -p "$RESULTS_DIR" "$FASTQ_TRIM_DIR" "$FASTQC_TRIM_DIR" 
mkdir -p "$ARRIBA_DIR" "$STAR_FUSION_DIR" "$REPORTS_DIR" "$FEATURE_COUNTS_DIR"

# Number of jobs to run in parallel, must be less than the number of samples
# Default to 1 if not provided
export PARALLEL_JOBS="${3:-8}"

# Specific for STAR align jobs, 1 jobs take 16 threads
export STAR_JOBS="${4:-1}"
export STAR_THREADS=16

# "=========================================================================="
# Reference and annotation directories
# "=========================================================================="
export REFERENCE_DIR="/mnt/f/projects/Reference"
# export STAR_INDEX="${REFERENCE_DIR}/Gencode/STAR_index"
export STAR_INDEX="${REFERENCE_DIR}/Gencode/STAR_index_hg38.v44"

# export REFERENCE="${REFERENCE_DIR}/Gencode/gencode.hg38.v36.primary_assembly.fa"
export REFERENCE="${REFERENCE_DIR}/Gencode/gencode.hg38.v44/GRCh38.primary_assembly.genome.fa"

# export ANNOTATION="${REFERENCE_DIR}/Gencode/gencode.v36.primary_assembly.annotation.gtf"
export ANNOTATION="${REFERENCE_DIR}/Gencode/gencode.hg38.v44/gencode.v44.primary_assembly.annotation.gtf"

export BAIT_INTERVALS="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2/hg38/xgen-exome-hyb-panel-v2-probes-hg38.interval_list"
export TARGET_INTERVALS="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2/hg38/xgen-exome-hyb-panel-v2-targets-hg38.interval_list"

export REF_FLAT="${REFERENCE_DIR}/Picard_QC/CollectRnaSeqMetrics/refFlat.txt"
export RIBO_INTERVALS="${REFERENCE_DIR}/Picard_QC/CollectRnaSeqMetrics/v44/GRCh38_gencode_v44_rRNA.interval_list"

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
echo "STAR align parallel jobs:     $STAR_JOBS (Default: 1)"
echo "STAR align job threads:       $STAR_THREADS (Default: 16)"
echo "Other steps parallel jobs:    $PARALLEL_JOBS (Default: 8)"
echo "======================================================================="
