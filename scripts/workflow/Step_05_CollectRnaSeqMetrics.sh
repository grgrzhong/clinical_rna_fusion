#!/bin/bash

basedir="/media/hkusarc/Data/STUMP_RNA_seq_ALL/Output_StarFusion/rerun" # Modify if needed
picard_jar="/media/hkusarc/Data/Software/picard/picard_3.3.0.jar"
ref_flat="/media/hkusarc/Data/Reference/Picard_QC/CollectRnaSeqMetrics/refFlat.txt"
ribo_intervals="/media/hkusarc/Data/Reference/Picard_QC/CollectRnaSeqMetrics/v44/GRCh38_gencode_v44_rRNA.interval_list"

# basedir="/media/hkusarc/Data/STUMP_RNA_seq_ALL/Output_StarFusion/rerun"  # Modify if needed
# picard_jar="/media/hkusarc/Data/Software/picard/picard_3.3.0.jar"
# REF_FLAT="/media/hkusarc/Data/Reference/Picard_QC/CollectRnaSeqMetrics/refFlat.txt"
# RIBO_INTERVALS="/media/hkusarc/Data/Reference/Picard_QC/CollectRnaSeqMetrics/v44/GRCh38_gencode_v44_rRNA.interval_list"

find "$basedir" -type f -name "Aligned.sortedByCoord.out.bam" | while read bam; do
    sample_dir=$(dirname "$bam")
    sample_name=$(basename "$sample_dir")
    output_metrics="$basedir/$sample_name/${sample_name}_RnaSeqMetrics.txt"
    output_dir=$(dirname "$output_metrics")

    mkdir -p "$output_dir" # Ensure the directory exists
    echo "Processing sample: $sample_name"

    java -jar "$picard_jar" CollectRnaSeqMetrics \
        I="$bam" \
        O="$output_metrics" \
        REF_FLAT="$REF_FLAT" \
        STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
        RIBOSOMAL_INTERVALS="$RIBO_INTERVALS"
done
