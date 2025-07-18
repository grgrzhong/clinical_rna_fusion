#!/bin/bash

#############################################################################
# Clinical RNA Fusion Analysis Workflow - STAR-Fusion Step
# This script processes RNA-Seq data using STAR-Fusion for fusion detection.
#############################################################################

# Create the output directory if it does not exist
mkdir -p "${STAR_FUSION_DIR}"

# Function to run STAR-Fusion for fusion detection 
star_fusion() {
    local sample="$1"

    ## Create output directory for STAR-Fusion results
    output_dir=${STAR_FUSION_DIR}/${sample}/
    mkdir -p "$output_dir"
    
    ## Temporary directory for STAR
    tmp_dir="/tmp/star_${sample}"
    rm -rf "$tmp_dir"

    ## Run STAR alignment for STAR-Fusion
    echo "$(date +"%F") $(date +"%T") - (${sample}) Running STAR alignment for STAR-Fusion ..."

    file1="${FASTQ_TRIM_DIR}/${sample}/${sample}_trimmed_R1.fastq.gz"
    file2="${FASTQ_TRIM_DIR}/${sample}/${sample}_trimmed_R2.fastq.gz"

    singularity exec \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind "${FASTQ_TRIM_DIR}:${FASTQ_TRIM_DIR}" \
        --bind "${STAR_FUSION_DIR}:${STAR_FUSION_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/star-fusion.sif" \
        STAR --genomeDir "${STAR_INDEX}" \
        --readFilesIn "$file1" "$file2" \
        --outReadsUnmapped None \
        --runThreadN "${STAR_THREADS}" \
        --twopassMode Basic \
        --readFilesCommand "gunzip -c" \
        --outSAMstrandField intronMotif \
        --outSAMunmapped Within \
        --chimSegmentMin 12 \
        --chimJunctionOverhangMin 8 \
        --chimOutJunctionFormat 1 \
        --alignSJDBoverhangMin 10 \
        --alignMatesGapMax 100000 \
        --alignIntronMax 100000 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --outSAMattrRGline ID:GRPundef SM:"$sample" \
        --chimMultimapScoreRange 3 \
        --chimScoreJunctionNonGTAG -4 \
        --chimMultimapNmax 20 \
        --chimNonchimScoreDropMin 10 \
        --peOverlapNbasesMin 12 \
        --peOverlapMMp 0.1 \
        --alignInsertionFlush Right \
        --alignSplicedMateMapLminOverLmate 0 \
        --alignSplicedMateMapLmin 30 \
        --outFileNamePrefix "${output_dir}" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outTmpDir "$tmp_dir" \
        >& "${output_dir}/star_star_align.log"

    # Run samtools with singularity exec
    echo "$(date +"%F") $(date +"%T") - (${sample}) Indexing BAM file with samtools ..."

    singularity exec \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind "${FASTQ_TRIM_DIR}:${FASTQ_TRIM_DIR}" \
        --bind "${STAR_FUSION_DIR}:${STAR_FUSION_DIR}" \
        "${CONTAINER_DIR}/star-fusion.sif" \
        samtools index "${output_dir}/Aligned.sortedByCoord.out.bam"

    # Run STAR-Fusion with singularity exec
    echo "$(date +"%F") $(date +"%T") - (${sample}) Running STAR Fusion ..."

    singularity exec \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind "${FASTQ_TRIM_DIR}:${FASTQ_TRIM_DIR}" \
        --bind "${STAR_FUSION_DIR}:${STAR_FUSION_DIR}" \
        "${CONTAINER_DIR}/star-fusion.sif" \
        STAR-Fusion \
        --genome_lib_dir "$CTAT_RESOURCE_LIB" \
        -J "$output_dir/Chimeric.out.junction" \
        --left_fq "$file1" \
        --right_fq "$file2" \
        --output_dir "${output_dir}" \
        --examine_coding_effect \
        --extract_fusion_reads \
        --FusionInspector inspect \
        >& "${output_dir}/star_detect_fusion.log"
}

# Export the function to be used by GNU parallel
export -f star_fusion

# Process samples in parallel
samples=$(find "${FASTQ_TRIM_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n")

echo "$samples" |
    parallel \
        --jobs "$STAR_JOBS" \
        star_fusion {}