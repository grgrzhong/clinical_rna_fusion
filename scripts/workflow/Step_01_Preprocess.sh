#!/bin/bash

## Function to process a single sample
preprocess() {
    local sample="$1"

    echo "$(date +"%F") $(date +"%T")" " - Processing sample = ${sample}"

    ## Create output directories
    mkdir -p "${FASTQ_TRIM_DIR}/${sample}"
    mkdir -p "${FASTQC_TRIM_DIR}/${sample}"

    ## Trim the fastq files using fastp
    singularity exec \
        --bind "${PROJECT_DIR}":"${PROJECT_DIR}" \
        --bind "${PRIMARY_SEQ_DIR}":"${PRIMARY_SEQ_DIR}" \
        --bind "${FASTQ_TRIM_DIR}":"${FASTQ_TRIM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/fastp.sif" \
        fastp \
        -i "${PRIMARY_SEQ_DIR}/${sample}_1.fastq.gz" \
        -I "${PRIMARY_SEQ_DIR}/${sample}_2.fastq.gz" \
        -o "${FASTQ_TRIM_DIR}/${sample}/${sample}_trimmed_R1.fastq.gz" \
        -O "${FASTQ_TRIM_DIR}/${sample}/${sample}_trimmed_R2.fastq.gz" \
        --detect_adapter_for_pe \
        --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --dont_eval_duplication \
        --trim_poly_g \
        --trim_poly_x \
        -j "${FASTQ_TRIM_DIR}/${sample}/${sample}.json" \
        -h "${FASTQ_TRIM_DIR}/${sample}/${sample}.html" \
        -w 2 \
        2>"${FASTQ_TRIM_DIR}/${sample}/${sample}.fastp.log"

    ## FastQC on trimmed fastq files
    echo "$(date +"%F") $(date +"%T")" " - Running FastQC on trimmed fastq files for sample ${sample}"
    singularity exec \
        --bind "${PROJECT_DIR}":"${PROJECT_DIR}" \
        --bind "${FASTQ_TRIM_DIR}":"${FASTQ_TRIM_DIR}" \
        --bind "${FASTQC_TRIM_DIR}":"${FASTQC_TRIM_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/fastqc.sif" \
        fastqc \
        "${FASTQ_TRIM_DIR}/${sample}/${sample}_trimmed_R1.fastq.gz" \
        "${FASTQ_TRIM_DIR}/${sample}/${sample}_trimmed_R2.fastq.gz" \
        -o "${FASTQC_TRIM_DIR}/${sample}" \
        --memory 4096 \
        -t 4 \
        >& "${FASTQC_TRIM_DIR}/${sample}/${sample}.fastqc.log"

    echo "$(date +"%F") $(date +"%T")" " - Completed preprocessing sample = ${sample}"
}

## Export the function so parallel can use it
export -f preprocess

# Get unique sample names from fastq files
samples=$(find "${PRIMARY_SEQ_DIR}" -name "*.fastq.gz" |
    sed 's/_[12]\.fastq\.gz$//' |
    sort -u |
    xargs -n1 basename)

# Process samples in parallel
echo "$samples" |
    parallel \
        --jobs "$JOBS" \
        --progress \
        preprocess {}
