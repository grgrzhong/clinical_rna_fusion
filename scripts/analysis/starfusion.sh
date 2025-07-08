#!/bin/bash
#SBATCH --job-name=RNA_STARFusion
#SBATCH --partition=amd
#SBATCH --time=48:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/home/zhonggr/projects/250224_DFSP_WES/slurm/%x_%j.out
#SBATCH --error=/home/zhonggr/projects/250224_DFSP_WES/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

##############################################################################
## Run STAR-Fusion using nextflow
##############################################################################
## Setup the working directory
# cd /home/zhonggr/projects/250224_DFSP_WES

# Global nextflow configuration
# export NXF_OFFLINE=true
# export NXF_HOME=$HOME/.nextflow
# export NXF_DISABLE_CHECK_LATEST=true
# export NXF_OPTS="-Xms512m -Xmx8g"
# export NXF_LOG_FILE="${PWD}/.nextflow.log"
# export NXF_WORK="${PWD}/work"

# echo ${NXF_WORK}
# rm -f NXF_WORK
# rm -f .nextflow.log*

# # # Run the Nextflow pipeline in local mode
# nextflow run workflows/rna_fusion.nf \
#     -profile local \
#     -resume

##############################################################################
## Run STAR-Fusion using Singularity
##############################################################################
# Load the Singularity module
# input_trimmed_dir=/mnt/m/RNA-seq/STUMP/Input-trimmed
input_trimmed_dir=/mnt/f/projects/250702_RNA_fusion/data
# ref_dir=/mnt/m/Reference
ref_dir=/mnt/f/projects/250702_RNA_fusion/reference
output_dir=/mnt/f/projects/250702_RNA_fusion/outputs/stump/STAR-Fusion
singularity_dir=/mnt/f/projects/250702_RNA_fusion/containers

# singularity shell \
#     --bind ${ref_dir}:${ref_dir} \
#     --bind ${input_trimmed_dir}:${input_trimmed_dir} \
#     --bind ${output_dir}:${output_dir} \
#     ${singularity_dir}/star-fusion.sif

STARINDEX=${ref_dir}/Gencode/STAR_index/
REFERENCE=${ref_dir}/Gencode/gencode.hg38.v36.primary_assembly.fa
ANNOTATION=${ref_dir}/Gencode/gencode.v36.primary_assembly.annotation.gtf
CTAT_RESOURCE_LIB=${ref_dir}/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/

# Check if input files exist
input_files=$(ls ${input_trimmed_dir}/*.fastq.gz 2>/dev/null | grep "R1")

# singularity exec -e -B `pwd` -B /path/to/ctat_genome_lib_build_dir \
#         star-fusion-v$version.simg \
#         STAR-Fusion \
#         --left_fq reads_1.fq.gz \
#         --right_fq reads_2.fq.gz \
#         --genome_lib_dir /path/to/ctat_genome_lib_build_dir \
#         -O StarFusionOut \
#         --FusionInspector validate \
#         --examine_coding_effect \
#         --denovo_reconstruct

# file=/mnt/m/RNA-seq/STUMP/Input-trimmed/S25_trimmed_R1.fastq.gz

for file in $(ls ${input_trimmed_dir}/*.fastq.gz | grep "R1"); do 

    echo $file; 
    
    FILENAME=$(basename $file | cut -d "_" -f 1); 
    
    echo $FILENAME; 
    
    OUTPUT=${output_dir}/${FILENAME}/
    mkdir -p $OUTPUT
    echo $OUTPUT;

    rm -rf /tmp/star_${FILENAME}/

    ## Run STAR alignment 
    singularity exec \
        --bind ${ref_dir}:${ref_dir} \
        --bind ${input_trimmed_dir}:${input_trimmed_dir} \
        --bind ${output_dir}:${output_dir} \
        --bind /tmp:/tmp \
        ${singularity_dir}/star-fusion.sif \
        STAR --genomeDir $STARINDEX \
            --readFilesIn $file ${file//R1/R2} \
            --outReadsUnmapped None \
            --runThreadN 8 \
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
            --outSAMattrRGline ID:GRPundef SM:$FILENAME \
            --chimMultimapScoreRange 3 \
            --chimScoreJunctionNonGTAG -4 \
            --chimMultimapNmax 20 \
            --chimNonchimScoreDropMin 10 \
            --peOverlapNbasesMin 12 \
            --peOverlapMMp 0.1 \
            --alignInsertionFlush Right \
            --alignSplicedMateMapLminOverLmate 0 \
            --alignSplicedMateMapLmin 30 \
            --outFileNamePrefix $OUTPUT \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            --outTmpDir /tmp/star_${FILENAME}/ \
            >& ${OUTPUT}/staralignment.log

    # bam file needs SM tag for CTAT mutation

    # Run samtools with singularity exec
    singularity exec \
        --bind ${ref_dir}:${ref_dir} \
        --bind ${input_trimmed_dir}:${input_trimmed_dir} \
        --bind ${output_dir}:${output_dir} \
        ${singularity_dir}/star-fusion.sif \
        samtools index ${OUTPUT}/Aligned.sortedByCoord.out.bam

    # Run STAR-Fusion with singularity exec
    singularity exec \
        --bind ${ref_dir}:${ref_dir} \
        --bind ${input_trimmed_dir}:${input_trimmed_dir} \
        --bind ${output_dir}:${output_dir} \
        ${singularity_dir}/star-fusion.sif \
        STAR-Fusion --genome_lib_dir $CTAT_RESOURCE_LIB \
                    -J $OUTPUT/Chimeric.out.junction \
                    --left_fq $file \
                    --right_fq ${file//R1/R2} \
                    --output_dir $OUTPUT \
                    --examine_coding_effect \
                    --extract_fusion_reads \
                    --FusionInspector inspect \
                    >& ${OUTPUT}/starfusion.log
done
