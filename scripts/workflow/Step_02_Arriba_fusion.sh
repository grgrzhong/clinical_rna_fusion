#!/bin/bash

#############################################################################
# Clinical RNA Fusion Analysis Workflow - Arriba Step
# This script processes RNA-Seq data using Arriba for fusion detection.
#############################################################################

# Loop through all fastq.gz files in the trimmed fastq directory
arriba_fusion(){
    local sample="$1"

    echo "$(date +"%F") $(date +"%T")" "- Processing sample = ${sample}"

    ## Create output directory for STAR-Fusion results
    output_dir=${ARRIBA_DIR}/${sample}/
    mkdir -p "$output_dir"
    
    ## Temporary directory for STAR
    tmp_dir="/tmp/star_${sample}"
    rm -rf "$tmp_dir"

    ## Run STAR alignment for Arriba
    echo "$(date +"%F") $(date +"%T")" "- Running STAR alignment for Arriba ..."
    
    file1="${FASTQ_TRIM_DIR}/${sample}/${sample}_trimmed_R1.fastq.gz"
    file2="${FASTQ_TRIM_DIR}/${sample}/${sample}_trimmed_R2.fastq.gz"
    
    singularity exec \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind "${FASTQ_TRIM_DIR}:${FASTQ_TRIM_DIR}" \
        --bind "${ARRIBA_DIR}:${ARRIBA_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/star-fusion.sif" \
        STAR \
        --runThreadN "${THREADS}" \
        --genomeDir "${STAR_INDEX}" \
        --readFilesIn "${file1}" "${file2}" \
        --readFilesCommand "gunzip -c" \
        --outFileNamePrefix "${output_dir}" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outBAMcompression 0 \
        --outFilterMultimapNmax 50 \
        --peOverlapNbasesMin 10 \
        --alignSplicedMateMapLminOverLmate 0.5 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimSegmentMin 10 \
        --chimOutType WithinBAM HardClip \
        --chimJunctionOverhangMin 10 \
        --chimScoreDropMax 30 \
        --chimScoreJunctionNonGTAG 0 \
        --chimScoreSeparation 1 \
        --chimSegmentReadGapMax 3 \
        --outTmpDir "${tmp_dir}" \
        --chimMultimapNmax 50 \
        >& "${output_dir}/arriba_star_align.log"

    ## Index the BAM file
    echo "$(date +"%F") $(date +"%T")" "- Indexing BAM file with samtools ..."
    rm -rf "${output_dir}/Aligned.sortedByCoord.out.bam.bai"
    
    singularity exec \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind "${FASTQ_TRIM_DIR}:${FASTQ_TRIM_DIR}" \
        --bind "${ARRIBA_DIR}:${ARRIBA_DIR}" \
        "${CONTAINER_DIR}/star-fusion.sif" \
        samtools index "${output_dir}/Aligned.sortedByCoord.out.bam"

    ## Run Arriba for fusion detection
    echo "$(date +"%F") $(date +"%T")" "- Running Arriba for fusion detection ..."
    
    singularity exec \
        --bind "${REFERENCE_DIR}":"${REFERENCE_DIR}" \
        --bind "${ARRIBA_DIR}":"${ARRIBA_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/arriba.sif" \
        arriba \
        -x "${output_dir}/Aligned.sortedByCoord.out.bam" \
        -o "${output_dir}/fusions.tsv" \
        -O "${output_dir}/fusions.discarded.tsv" \
        -a "${REFERENCE}" \
        -g "${ANNOTATION}" \
        -b "${BLACKLIST}" \
        -k "${KNOWN_FUSION}" \
        -t "${KNOWN_FUSION}" \
        -p "${PROTEIN_DOMAINS}" \
        >& "${output_dir}/arriba_detect_fusion.log"

    ## Draw Arriba fusions
    echo "$(date +"%F") $(date +"%T")" "- Drawing Arriba fusions ..."
    
    singularity exec \
        --bind "${REFERENCE_DIR}":"${REFERENCE_DIR}" \
        --bind "${ARRIBA_DIR}":"${ARRIBA_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/arriba.sif" \
        draw_fusions.R \
        --fusions="${output_dir}/fusions.tsv" \
        --output="${output_dir}/fusions.pdf" \
        --annotation="${ANNOTATION}" \
        --cytobands="${CYTOBANDS}" \
        --proteinDomains="${PROTEIN_DOMAINS}" \
        --minConfidenceForCircosPlot=none \
        --mergeDomainsOverlappingBy=0.5 \
        >& "${output_dir}/arriba_draw_fusions.log"

    ## Keep only files start with "fusions" and log files in the output directory
    echo "$(date +"%F") $(date +"%T")" "- Keeping only relevant files in output directory for Arriba ..."
    find "$output_dir" -type f ! -name "fusions*" ! -name "*.log" -exec rm -f {} \;

}

# Export the function for parallel execution
export -f arriba_fusion

# Process samples in parallel
samples=$(find "${FASTQ_TRIM_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n")

echo "$samples" |
    parallel \
        --jobs "$JOBS" \
        --progress \
        arriba_fusion {}