#!/bin/bash

## Funtion to generate the expression matrix
clean_remove_duplicates() {
    local sample="$1"

    echo "$(date +"%F") $(date +"%T") - Processing sample = ${sample}"

    ## Sort the BAM file
    echo "$(date +"%F") $(date +"%T") - Sorting BAM file ..."
    input_bam="${STAR_FUSION_DIR}/${sample}/Aligned.sortedByCoord.out.bam"

    ## Check if input BAM file exists
    if [[ ! -f "$input_bam" ]]; then
        echo "Warning: Input BAM file not found: $input_bam"
        return 1
    fi

    ## Create output directories if they do not exist
    mkdir -p "${FEATURE_COUNTS_DIR}/${sample}"

    clean_bam="${FEATURE_COUNTS_DIR}/${sample}/${sample}.clean.bam"
    dedup_bam="${FEATURE_COUNTS_DIR}/${sample}/${sample}.dedup.bam"
    metrics_file="${FEATURE_COUNTS_DIR}/${sample}/${sample}.dedup.metrics.txt"

    ## Filter for properly paired, unique, high-quality reads
    echo "$(date +"%F") $(date +"%T") - Filtering BAM file for properly paired, unique, high-quality reads ..."
    singularity exec \
        --bind "${STAR_FUSION_DIR}:${STAR_FUSION_DIR}" \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind "${FEATURE_COUNTS_DIR}:${FEATURE_COUNTS_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/samtools.sif" \
        bash -c "samtools view -h -q 20 -f 0x2 -F 0x904 '${input_bam}' | awk '\$0 ~ /^@/ || \$0 ~ /NH:i:1/' | samtools sort -o '${clean_bam}'"

    ## Mark and remove PCR duplicates
    echo "$(date +"%F") $(date +"%T") - Marking and removing PCR duplicates for sample = ${sample}"
    singularity exec \
        --bind "${STAR_FUSION_DIR}:${STAR_FUSION_DIR}" \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind "${FEATURE_COUNTS_DIR}:${FEATURE_COUNTS_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/picard.sif" \
        picard "-Xmx4g" \
        MarkDuplicates \
        I="${clean_bam}" \
        O="${dedup_bam}" \
        M="${metrics_file}" \
        REMOVE_DUPLICATES=true

    ## Index the deduplicated BAM file
    echo "$(date +"%F") $(date +"%T") - Indexing deduplicated BAM file ..."
    singularity exec \
        --bind "${STAR_FUSION_DIR}:${STAR_FUSION_DIR}" \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind "${FEATURE_COUNTS_DIR}:${FEATURE_COUNTS_DIR}" \
        "${CONTAINER_DIR}/samtools.sif" \
        samtools index "${dedup_bam}"

    ## Clean up intermediate files
    echo "$(date +"%F") $(date +"%T") - Cleaning up intermediate files ..."
    rm -rf "${clean_bam}"
}

# Export the function for parallel execution
export -f clean_remove_duplicates

# Process samples in parallel
samples=$(find "${STAR_FUSION_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n")

echo "$samples" |
    parallel \
        --jobs "$PARALLEL_JOBS" \
        clean_remove_duplicates {}

# Generate feature counts for expression matrix
num_samples=$(echo "$samples" | wc -l)
echo "$(date +"%F") $(date +"%T") - Generating featureCounts expression matrix for ${num_samples} samples"

bam_files=$(find "${FEATURE_COUNTS_DIR}" -type f -name "*dedup.bam" ! -name "*dedup.bam.bai" | tr '\n' ' ')

echo "$bam_files"

feature_counts_file="${FEATURE_COUNTS_DIR}/merged.featureCounts.txt"

singularity exec \
    --bind "${STAR_FUSION_DIR}:${STAR_FUSION_DIR}" \
    --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
    --bind "${FEATURE_COUNTS_DIR}:${FEATURE_COUNTS_DIR}" \
    --bind /tmp:/tmp \
    "${CONTAINER_DIR}/subread.sif" \
    featureCounts -p -B -C \
    -a "$ANNOTATION" \
    -o "${feature_counts_file}" \
    -s 2 \
    -T 16 \
    -t exon \
    ${bam_files} #!! Should not be quotes when using singularity
