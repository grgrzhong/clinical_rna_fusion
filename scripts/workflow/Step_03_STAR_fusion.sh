#!/bin/bash


## Loop through all fastq.gz files in the trimmed fastq directory
for file in $(ls ${FASTQ_TRIM_DIR}/*.fastq.gz | grep "R1"); do

    echo "$(date +"%F") $(date +"%T")" "Processing sample = ${file}"

    sample=$(basename "$file" | cut -d "_" -f 1)

    ## Create output directory for STAR-Fusion results
    output_dir=${STAR_FUSION_DIR}/${sample}/
    mkdir -p "$output_dir"
    
    ## Temporary directory for STAR
    tmp_dir="/tmp/star_${sample}"
    rm -rf "$tmp_dir"

    ## Run STAR alignment for STAR-Fusion
    echo "$(date +"%F") $(date +"%T")" "- Running STAR alignment for STAR-Fusion ..."
    
    singularity exec \
        --bind ${REFERENCE_DIR}:${REFERENCE_DIR} \
        --bind ${FASTQ_TRIM_DIR}:${FASTQ_TRIM_DIR} \
        --bind ${STAR_FUSION_DIR}:${STAR_FUSION_DIR} \
        --bind /tmp:/tmp \
        ${CONTAINER_DIR}/star-fusion.v1.15.0.simg \
        STAR --genomeDir "${STAR_INDEX}" \
        --readFilesIn "$file" "${file//R1/R2}" \
        --outReadsUnmapped None \
        --runThreadN "${THREADS}" \
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
    echo "$(date +"%F") $(date +"%T")" "- Indexing BAM file with samtools ..."

    singularity exec \
        --bind ${REFERENCE_DIR}:${REFERENCE_DIR} \
        --bind ${FASTQ_TRIM_DIR}:${FASTQ_TRIM_DIR} \
        --bind ${STAR_FUSION_DIR}:${STAR_FUSION_DIR} \
        "${CONTAINER_DIR}/star-fusion.v1.15.0.simg" \
        samtools index "${output_dir}/Aligned.sortedByCoord.out.bam"

    # Run STAR-Fusion with singularity exec
    echo "$(date +"%F") $(date +"%T")" "- Running STAR Fusion ..."

    singularity exec \
        --bind ${REFERENCE_DIR}:${REFERENCE_DIR} \
        --bind ${FASTQ_TRIM_DIR}:${FASTQ_TRIM_DIR} \
        --bind ${STAR_FUSION_DIR}:${STAR_FUSION_DIR} \
        "${CONTAINER_DIR}/star-fusion.v1.15.0.simg" \
        STAR-Fusion \
        --genome_lib_dir "$CTAT_RESOURCE_LIB" \
        -J "$output_dir/Chimeric.out.junction" \
        --left_fq "$file" \
        --right_fq "${file//R1/R2}" \
        --output_dir "${output_dir}" \
        --examine_coding_effect \
        --extract_fusion_reads \
        --FusionInspector inspect \
        >& "${output_dir}/star_detect_fusion.log"
done
