#!/bin/bash

## Loop thorough all fastq.gz files in the primary sequence directory
for file in $(find "${PRIMARY_SEQ_DIR}" -name "*.fastq.gz" | rev | cut -d '_' -f 2- | rev | uniq); do

    echo "$(date +"%F") $(date +"%T")" "Processing sample = ${file}"

    sample=$(basename "$file")

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
        -o "${FASTQ_TRIM_DIR}/${sample}_trimmed_R1.fastq.gz" \
        -O "${FASTQ_TRIM_DIR}/${sample}_trimmed_R2.fastq.gz" \
        --detect_adapter_for_pe \
        --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --dont_eval_duplication \
        --trim_poly_g \
        --trim_poly_x \
        -j "${FASTQ_TRIM_DIR}/${sample}.json" \
        -h "${FASTQ_TRIM_DIR}/${sample}.html" \
        -w 8 \
        2>"${FASTQ_TRIM_DIR}/${sample}.fastp.log"

done

## FastQC on trimmed fastq files
echo "$(date +"%F") $(date +"%T")" "Running FastQC on trimmed fastq files"

singularity exec \
    --bind "${PROJECT_DIR}":"${PROJECT_DIR}" \
    --bind "${FASTQ_TRIM_DIR}":"${FASTQ_TRIM_DIR}" \
    --bind "${FASTQC_TRIM_DIR}":"${FASTQC_TRIM_DIR}" \
    --bind /tmp:/tmp \
    "${CONTAINER_DIR}/fastqc.sif" \
    fastqc \
    "${FASTQ_TRIM_DIR}"/*.fastq.gz \
    -o "${FASTQC_TRIM_DIR}" \
    --memory 8192 \
    -t 8 \
    2>"${FASTQC_TRIM_DIR}/fastqc.log"
