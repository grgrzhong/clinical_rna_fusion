#!/bin/bash

##### prerequisite #####
# conda activate fusion_report_env
# Python version requirement 3.7

##### test #####
# Run the fusion-report - one sample test
# fusion_report run "RNA195" /media/hkusarc/Data/clinical_rnaseq_june11/fusion_report /media/hkusarc/Data/Reference/fusion_report \
#   --arriba /media/hkusarc/Data/clinical_rnaseq_june11/Output-exome-Arriba/RNA195/fusions.tsv \
#   --starfusion /media/hkusarc/Data/clinical_rnaseq_june11/Output-exome-STARFusion/RNA195/star-fusion.fusion_predictions.tsv \
#   --no-cosmic \
#   --allow-multiple-gene-symbols

##### fusion_report #####
# Run the fusion-report - all sample
# Base directories
base_dir="/mnt/f/projects/250702_RNA_fusion"
arriba_dir="$base_dir/outputs/Arriba"
starfusion_dir="$base_dir/outputs/StarFusion"
reference_dir="/mnt/f/projects/Reference/fusion_report"
container_dir="/mnt/f/projects/250702_RNA_fusion/containers"

# Iterate over all sample directories in STARFusion directory
for sample_dir in "$starfusion_dir"/*; do
    if [ -d "$sample_dir" ]; then
        # Extract sample name dynamically
        sample_name=$(basename "$sample_dir")

        # Define input files
        arriba_file="$arriba_dir/$sample_name/fusions.tsv"
        starfusion_file="$starfusion_dir/$sample_name/star-fusion.fusion_predictions.tsv"

        # Extract sample names from file paths
        arriba_sample_name=$(basename "$(dirname "$arriba_file")")
        starfusion_sample_name=$(basename "$(dirname "$starfusion_file")")

        # Validate that the sample name matches between ARRIBA and STARFusion files
        if [[ "$arriba_sample_name" == "$sample_name" && "$starfusion_sample_name" == "$sample_name" ]]; then
            # Define output directory for the sample
            sample_output_dir="$starfusion_dir/$sample_name/fusion_report"
            mkdir -p "$sample_output_dir"

            # Run fusion_report for the sample
            singularity exec \
                --bind ${base_dir}:${base_dir} \
                --bind ${arriba_dir}:${arriba_dir} \
                --bind ${starfusion_dir}:${starfusion_dir} \
                --bind ${reference_dir}:${reference_dir} \
                --bind /tmp:/tmp \
                ${container_dir}/fusion_report.sif fusion_report run "$sample_name" "$sample_output_dir" "$reference_dir" \
                --arriba "$arriba_file" \
                --starfusion "$starfusion_file" \
                --no-cosmic \
                --allow-multiple-gene-symbols
        else
            echo "Warning: Missing input files for sample $sample_name. Skipping..."
        fi
    fi
done
