#!/bin/bash
#SBATCH --job-name=clinical_rna_fusion
#SBATCH --partition=amd
#SBATCH --time=6:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=3G
#SBATCH --output=/lustre1/g/path_my/pipeline/clinical_rna_fusion/slurm/%x_%j.out
#SBATCH --error=/lustre1/g/path_my/pipeline/clinical_rna_fusion/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

## Reference: https://github.com/deweylab/RSEM
## Create the environment if not already created
source $(conda info --base)/etc/profile.d/conda.sh
# conda create -n rsem -c bioconda -c conda-forge fastp=0.23.4 star=2.7.10b rsem=1.3.3 samtools=1.20 parallel
conda activate rsem

## Create the output directory if it does not exist
export REFERENCE_DIR=/lustre1/g/path_my/Reference
export REFERENCE="${REFERENCE_DIR}/Gencode/gencode.hg38.v44/GRCh38.primary_assembly.genome.fa"
export ANNOTATION="${REFERENCE_DIR}/Gencode/gencode.hg38.v44/gencode.v44.primary_assembly.annotation.gtf"
# export ANNOTATION="${REFERENCE_DIR}/Gencode/gencode.hg38.v44/gencode.v44.protein_coding.annotation.gtf"
export OUTPUT_DIR=/lustre1/g/path_my/pipeline/clinical_rna_fusion/data/DFSP
export FASTQ_TRIM_DIR="${OUTPUT_DIR}/Input-trimmed"

## Create a map of ENSEMBL gene IDs to transcript IDs and gene names
awk -F'\t' '$3=="transcript" {
    # Extract gene_id, transcript_id, and gene_name from the 9th field
    gene_id = ""; transcript_id = ""; gene_name = "";
    split($9, attrs, ";");
    for (i in attrs) {
        if (match(attrs[i], /gene_id "([^"]+)"/, arr)) {
            gene_id = arr[1];
            # Remove version number (everything after the dot)
            gsub(/\.[0-9]+$/, "", gene_id);
        }
        if (match(attrs[i], /transcript_id "([^"]+)"/, arr)) {
            transcript_id = arr[1];
            # Remove version number (everything after the dot)
            gsub(/\.[0-9]+$/, "", transcript_id);
        }
        if (match(attrs[i], /gene_name "([^"]+)"/, arr)) {
            gene_name = arr[1];
        }
    }
    if (gene_id != "" && transcript_id != "" && gene_name != "") {
        print gene_id "\t" transcript_id "\t" gene_name;
    }
}' "${ANNOTATION}" > "${REFERENCE_DIR}/Gencode/gencode.hg38.v44/GRCh38.primary_assembly.genecode.v44.id_map.txt"

# export STAR_PATH=~/miniforge3/envs/rnaseq/bin/STAR
# export STAR_FUSION_DIR="${OUTPUT_DIR}/Output"
# export STAR_INDEX="${REFERENCE_DIR}/Gencode/STAR_index_hg38.v44"
export THREADS=16
export JOBS=2

export RSEM_DIR="${OUTPUT_DIR}/RSEM" && mkdir -p "$RSEM_DIR"
export RSEM_REF_DIR="${REFERENCE_DIR}/RSEM"
export RSEM_REF_NAME="GRCh38_gencode_v44"

## "======================================================================="
## Generate the RSEM reference
## "======================================================================="
# if [[ ! -f "${RSEM_REF_DIR}/${RSEM_REF_NAME}.grp" ]]; then
    echo "$(date +"%F") $(date +"%T") Preparing RSEM reference ..."
    mkdir -p "$RSEM_REF_DIR"

    # Build RSEM reference
    rsem-prepare-reference \
        --gtf "${ANNOTATION}" \
        --star \
        -p 16 \
        "${REFERENCE}" \
        "${RSEM_REF_DIR}/${RSEM_REF_NAME}" \
        >& "${RSEM_REF_DIR}/rsem_prepare_reference.log"

# rsem-prepare-reference \
# --gtf "/lustre1/g/path_my/Reference/Gencode/gencode.hg38.v44/gencode.v44.primary_assembly.annotation.gtf" \
# -p 8 \
# "/lustre1/g/path_my/Reference/Gencode/gencode.hg38.v44/GRCh38.primary_assembly.genome.fa" \
# "/lustre1/g/path_my/Reference/RSEM/GRCh38_gencode_v44"
# fi

## "======================================================================="
## STAR alignment for RSEM quantification
## "======================================================================="
star_alignment() {
    
    local sample="$1"
    # sample=DFSP-010-T-P1

    ## Create output directory for RSEM results
    output_dir="${RSEM_DIR}/${sample}"
    mkdir -p "$output_dir"
    
    ## Input files
    file1="${FASTQ_TRIM_DIR}/${sample}/${sample}_trimmed_R1.fastq.gz"
    file2="${FASTQ_TRIM_DIR}/${sample}/${sample}_trimmed_R2.fastq.gz"
    
    ## Check if input files exist
    if [[ ! -f "$file1" || ! -f "$file2" ]]; then
        echo "$(date +"%F") $(date +"%T") - ERROR: Input files not found for sample ${sample}"
        return 1
    fi
        
    ## Run STAR alignment for RSEM (transcriptome alignment)
    echo "$(date +"%F") $(date +"%T") - (${sample}) Running STAR alignment for RSEM ..."
    
    STAR \
        --genomeDir "${STAR_INDEX}" \
        --readFilesIn "$file1" "$file2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${output_dir}/${sample}.STAR." \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --sjdbScore 1 \
        --runThreadN "${THREADS}" \
        --quantMode TranscriptomeSAM \
        --outSAMattrRGline ID:"$sample" SM:"$sample" \
        >& "${output_dir}/star_alignment.log"
    
    ## sort genome bam
    samtools sort -@ "${THREADS}" -O BAM --write-index \
    -o "${output_dir}/${sample}.STAR.Aligned.toGenome.sort.bam" \
    "${output_dir}/${sample}.STAR.Aligned.out.bam"

    ## Remove unsorted genome bam file
    if [[ -e "${output_dir}/${sample}.STAR.Aligned.out.bam" ]]; then
        rm "${output_dir}/${sample}.STAR.Aligned.out.bam"
    fi

}

export -f star_alignment


## "======================================================================="
## RSEM quantification using the STAR-aligned transcriptome BAM
## "======================================================================="
rsem_quantification() {

    local sample="$1"

    ## Create output directory for RSEM results
    output_dir="${RSEM_DIR}/${sample}"
    mkdir -p "$output_dir"
    
    echo "$(date +"%F") $(date +"%T") - (${sample}) Running RSEM quantification ..."
    
    rsem-calculate-expression \
        --alignments \
        --no-bam-output \
        --paired-end \
        --forward-prob 0.5 \
        -p "${THREADS}" \
        "${output_dir}/${sample}.STAR.Aligned.toTranscriptome.out.bam" \
        "${RSEM_REF_DIR}/${RSEM_REF_NAME}" \
        "${output_dir}/${sample}.rsem" \
        >& "${output_dir}/rsem_quantification.log"
}

export -f rsem_quantification

## "======================================================================="
## Process samples in parallel
## "======================================================================="
export samples=$(find "${FASTQ_TRIM_DIR}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n")

# echo "$samples" | parallel --jobs "$JOBS" star_alignment {}
# echo "$samples" | parallel --jobs "$JOBS" rsem_quantification {}