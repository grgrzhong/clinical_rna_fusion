#!/bin/bash

RESULTS_DIR=/mnt/f/projects/clinical_rna_fusion
export CONTAINER_DIR="${PROJECT_DIR}/containers"
export MODULE_DIR="${PROJECT_DIR}/scripts/modules"
export RESULTS_DIR="/mnt/f/projects/clinical_rna_fusion/results/Exome"
export FASTQ_TRIM_DIR="${RESULTS_DIR}/Input-trimmed"
export FASTQC_TRIM_DIR="${RESULTS_DIR}/FastQC-trimmed"
export ARRIBA_DIR="${RESULTS_DIR}/Output"
export STAR_FUSION_DIR="${RESULTS_DIR}/Output"

rnafusion_pipeline.sh ./data/Raw ./results/Exome 8 1