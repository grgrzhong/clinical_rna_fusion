#!/bin/bash
exec > >(tee -i script.log)  # Log stdout
exec 2> >(tee -i script.err >&2)  # Log stderr

STARINDEX=/media/hkusarc/Data/Reference/Gencode/STAR_index_2024
CTAT_RESOURCE_LIB=/media/hkusarc/Data/Reference/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir
# INPUT=/media/hkusarc/Data/clinical_rnaseq_june11/Input-trimmed
INPUT=/media/hkusarc/Data/Input-trimmed

export PATH="/media/hkusarc/Data/Software/arriba_v2.4.0:/media/hkusarc/Data/Software/STAR-2.7.11b/bin/Linux_x86_64:$PATH"
export TRINITY_HOME=/home/hkusarc/miniconda3/envs/starfusion/opt/trinity-2.5.1

for file in $(ls "${INPUT}"/*.fastq.gz | grep "R1"); do 
    echo "Processing $file"
    FILENAME=$(basename "$file" | cut -d "_" -f 1)
    echo "Sample: $FILENAME"

    OUTPUT=$(dirname "$file")/../Output-exome-STARFusion/"$FILENAME"/
    mkdir -p "$OUTPUT"

    STAR --genomeDir "$STARINDEX" \
        --readFilesIn "$file" "${file//R1/R2}" \
        --outReadsUnmapped None \
        --runThreadN 16 \
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
        --outSAMattrRGline ID:GRPundef SM:"$FILENAME" \
        --chimMultimapScoreRange 3 \
        --chimScoreJunctionNonGTAG -4 \
        --chimMultimapNmax 20 \
        --chimNonchimScoreDropMin 10 \
        --peOverlapNbasesMin 12 \
        --peOverlapMMp 0.1 \
        --alignInsertionFlush Right \
        --alignSplicedMateMapLminOverLmate 0 \
        --alignSplicedMateMapLmin 30 \
        --outFileNamePrefix "$OUTPUT" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        >& "${OUTPUT}/staralignment.log"

    if [ -f "${OUTPUT}/Aligned.sortedByCoord.out.bam" ]; then
        samtools index "${OUTPUT}/Aligned.sortedByCoord.out.bam"
    else
        echo "Error: BAM file not found for $FILENAME"
        continue
    fi

    STAR-Fusion --genome_lib_dir "$CTAT_RESOURCE_LIB" \
                -J "${OUTPUT}/Chimeric.out.junction" \
                --left_fq "$file" \
                --right_fq "${file//R1/R2}" \
                --output_dir "$OUTPUT" \
                --examine_coding_effect \
                --extract_fusion_reads \
                --FusionInspector inspect \
                # --FusionInspector-opts "--min_junction_reads 1 --min_spanning_frags_only 0" \
                >& "${OUTPUT}/starfusion.log"
done