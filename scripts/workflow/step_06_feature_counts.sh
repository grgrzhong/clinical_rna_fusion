#!/bin/bash

## Funtion to generate the expression matrix
feature_counts() {
    local sample="$1"

    echo "$(date +"%F") $(date +"%T") - Processing sample = ${sample}"

    ## Sort the BAM file
    echo "$(date +"%F") $(date +"%T") - Sorting BAM file for sample = ${sample}"
    input_bam="${STAR_FUSION_DIR}/${sample}/Aligned.sortedByCoord.out.bam"
    
    ## Check if input BAM file exists
    if [[ ! -f "$input_bam" ]]; then
        echo "Warning: Input BAM file not found: $input_bam"
        return 1
    fi
    
    ## Generate the expression matrix using featureCounts
    echo "$(date +"%F") $(date +"%T") - Running featureCounts for sample = ${sample}"
    feature_count_file="${FEATURE_COUNTS_DIR}/${sample}/${sample}.featureCounts.txt"

    singularity exec \
        --bind "${STAR_FUSION_DIR}:${STAR_FUSION_DIR}" \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind "${FEATURE_COUNTS_DIR}:${FEATURE_COUNTS_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/subread.sif" \
        featureCounts \
        -p \
        -T 4 \
        -s 2 \
        -a "$ANNOTATION" \
        -o "${feature_count_file}" \
        "${input_bam}" \
        >& "${FEATURE_COUNTS_DIR}/${sample}/${sample}.featureCounts.log"
}

# Export the function for parallel execution
export -f feature_counts

# Process samples in parallel
samples=$(find "${STAR_FUSION_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n")

echo "$samples" |
    parallel \
        --jobs "$JOBS" \
        --progress \
        feature_counts {}