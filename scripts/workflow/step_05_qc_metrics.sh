#!/bin/bash

#############################################################################
# Clinical RNA Fusion Analysis Workflow - QC Metrics Step
# This script collects quality control metrics for RNA-Seq data using Picard
############################################################################

## Collect HS metrics for each sample in STAR_FUSION_DIR
generate_qc_metrics() {
    
    local sample="$1"

    input_bam="${STAR_FUSION_DIR}/${sample}/Aligned.sortedByCoord.out.bam"
    
    ## Check if input BAM file exists
    if [[ ! -f "$input_bam" ]]; then
        echo "Warning: Input BAM file not found: $input_bam"
        return 1
    fi
    
    # CollectHsMetrics information
    echo "$(date +"%F") $(date +"%T") - (${sample}) Running CollectHsMetrics ..."
    
    hs_metrics_file="${STAR_FUSION_DIR}/${sample}/${sample}_hs_metrics.txt"

    singularity exec \
        --bind "${STAR_FUSION_DIR}:${STAR_FUSION_DIR}" \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/picard.sif" \
        picard "-Xmx4g" \
        CollectHsMetrics \
        -I "${input_bam}" \
        -O "${hs_metrics_file}" \
        -R "${REFERENCE}" \
        -BI "${BAIT_INTERVALS}" \
        -TI "${TARGET_INTERVALS}" \
        >& "${STAR_FUSION_DIR}/${sample}/${sample}_hs_metrics.log"

    # CollectRnaSeqMetrics information
    echo "$(date +"%F") $(date +"%T") - (${sample}) Running CollectRnaSeqMetrics ..."

    rna_metrics_file="${STAR_FUSION_DIR}/${sample}/${sample}_rna_metrics.txt"

    singularity exec \
        --bind "${STAR_FUSION_DIR}:${STAR_FUSION_DIR}" \
        --bind "${REFERENCE_DIR}:${REFERENCE_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/picard.sif" \
        picard "-Xmx4g" \
        CollectRnaSeqMetrics \
        --INPUT "${input_bam}" \
        --OUTPUT "${rna_metrics_file}" \
        --REF_FLAT "$REF_FLAT" \
        --STRAND_SPECIFICITY SECOND_READ_TRANSCRIPTION_STRAND \
        --RIBOSOMAL_INTERVALS "$RIBO_INTERVALS" \
        >& "${STAR_FUSION_DIR}/${sample}/${sample}_rna_metrics.log"

    ## Generate MultiQC report for 
    echo "$(date +"%F") $(date +"%T") - (${sample}) Generating MultiQC report ..."
    
    multiqc_dir="${STAR_FUSION_DIR}/${sample}/multiqc"

    mkdir -p "$multiqc_dir"

    multiqc_filename="${sample}_QC_Report"

    singularity exec \
        --bind "${STAR_FUSION_DIR}:${STAR_FUSION_DIR}" \
        --bind /tmp:/tmp \
        "${CONTAINER_DIR}/multiqc.sif" \
        multiqc \
        --outdir "${multiqc_dir}" \
        --filename "${multiqc_filename}" \
        --title "QC Report - Sample ${sample}" \
        --comment "Quality control metrics for sample ${sample}" \
        --force \
        "${STAR_FUSION_DIR}/${sample}" \
        >& "${multiqc_dir}/${sample}_multiqc.log"
    
    ## Export the metrics to CSV
    echo "$(date +"%F") $(date +"%T") - (${sample}) Exporting metrics to CSV ..."

    singularity exec \
        --bind "${STAR_FUSION_DIR}:${STAR_FUSION_DIR}" \
        --bind "${multiqc_dir}:${multiqc_dir}" \
        "${CONTAINER_DIR}/r.sif" \
        Rscript "${MODULE_DIR}/combine_picard_hs_rna_metrics.R" \
        --hs_metrics "${multiqc_dir}/${multiqc_filename}_data/multiqc_picard_HsMetrics.txt" \
        --rna_metrics "${multiqc_dir}/${multiqc_filename}_data/multiqc_picard_RnaSeqMetrics.txt" \
        --sample_name "${sample}" \
        --output_csv "${STAR_FUSION_DIR}/${sample}/${sample}_RnaSeqMetrics.csv"
}

# Export the function for parallel execution
export -f generate_qc_metrics

# Process samples in parallel
samples=$(find "${STAR_FUSION_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n")

echo "$samples" |
    parallel \
        --jobs "$PARALLEL_JOBS" \
        generate_qc_metrics {}