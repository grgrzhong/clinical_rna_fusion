Clinical RNA Fusion Analysis Pipeline

![Pipeline Status](https://img.shields.io/badge/pipeline-active-brightgreen.svg)
![Platform](https://img.shields.io/badge/platform-Linux-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

A simple pipeline for detecting RNA fusions in clinical samples using Arriba and STAR-Fusion software. This pipeline processes paired-end RNA-seq FASTQ files through quality control, alignment, fusion detection, and generates detailed reports suitable for clinical applications.

- [Requirements](#requirements)
- [Pipeline details](#pipeline-details)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Resource](#resource)
- [Quick start](#quick-start)
- [Troubleshooting](#troubleshooting)

## Requirements

Install conda and essential software
```bash
## Download and install Miniforge (includes conda and mamba)
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh
## Follow installation prompts and restart terminal
source ~/.bashrc

## Verify installation
conda --version
mamba --version

## Create the base environemnt with apptainer and gnuparallel
conda env create -f ./conf/base.yml
```

## Pipeline details

| Step | Script | Tools Used | Description |
|------|--------|------------|-------------|
| **1** | `step_01_preprocess.sh` | fastp, FastQC | Quality control and adapter trimming |
| **2** | `step_02_arriba_fusion.sh` | STAR, Arriba | STAR alignment and Arriba fusion detection |
| **3** | `step_03_star_fusion.sh` | STAR-Fusion | Dedicated STAR-Fusion analysis |
| **4** | `step_04_fusion_report.sh` | Custom R scripts | Generate comprehensive fusion reports |
| **5** | `step_05_qc_metrics.sh` | Picard, MultiQC | Collect QC metrics and generate reports |

## Inputs
Run the pipeline with Key parameters:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `INPUT_DIR` | Raw FASTQ files directory | `${PROJECT_DIR}/data/test_data/Raw` |
| `OUTPUT_DIR` | Results output directory | `${PROJECT_DIR}/results` |
| `PARALLEL_JOBS` | Parallel sample processing | `1` |
| `STAR_JOBS` | STAR alignment jobs (16 threads each) | `1` |
| `REFERENCE_DIR` | Reference files directory | `/mnt/f/Reference` |
| `CONTAINER_DIR` | Singularity containers directory | `${PROJECT_DIR}/containers` |

The pipeline expects paired-end FASTQ files following this pattern:
```bash
{sample_id}_1.fastq.gz  # Forward reads (R1)
{sample_id}_2.fastq.gz  # Reverse reads (R2)
```
## Outputs

The pipeline generates a comprehensive output structure:

```
output_directory/
├── Input-trimmed/                    # Trimmed FASTQ files
│   └── {sample}/
│       ├── {sample}_trimmed_R1.fastq.gz
│       ├── {sample}_trimmed_R2.fastq.gz
│       ├── {sample}.json             # fastp report
│       └── {sample}.html             # fastp HTML report
├── FastQC-trimmed/                   # Quality control reports
│   └── {sample}/
│       ├── {sample}_trimmed_R1_fastqc.html
│       ├── {sample}_trimmed_R2_fastqc.html
│       └── {sample}.fastqc.log
├── Output/                           # Main analysis results
│   └── {sample}/
│       ├── Aligned.sortedByCoord.out.bam      # STAR alignment
│       ├── Chimeric.out.junction              # Chimeric reads
│       ├── fusions.tsv                        # Arriba results
│       ├── fusions.discarded.tsv              # Filtered fusions
│       ├── star-fusion.fusion_predictions.tsv # STAR-Fusion results
│       ├── {sample}_RnaSeqMetrics.txt         # Picard RNA metrics
│       ├── {sample}_HsMetrics.txt             # Hybrid selection metrics
│       └── fusion_report/                     # Combined reports
│           ├── fusion_report.xlsx
│           └── fusion_report.json
├── Reports/                          # Summary reports
│   ├── multiqc_report.html          # MultiQC summary
│   └── cohort_qc_metrics.csv        # Cohort-level metrics
└── Feature-counts/                   # Gene expression counts
    └── {sample}/
        └── {sample}_featurecounts.txt
```

## Resource

> ⚠️ **Important**: Ensure `STAR_JOBS × 16 ≤ total_CPU_cores` to prevent system overload.

For optimal performance:
- **Single sample**: `PARALLEL_JOBS=1, STAR_JOBS=1` (16 cores)
- **Multiple samples (32 cores)**: `PARALLEL_JOBS=4, STAR_JOBS=2` 
- **High-throughput (64+ cores)**: `PARALLEL_JOBS=8, STAR_JOBS=4`

## Quick start

```bash
# Basic usage with default test data
./run_pipeline_local.sh

# Custom parameters
./rnafusion_pipeline.sh \
    /path/to/fastq/files \
    /path/to/output \
    4 \
    2 \
    /path/to/reference \
    /path/to/containers
```

```bash
# Submit to SLURM scheduler
./run_pipeline_hpc.sh

# Or submit directly
sbatch rnafusion_pipeline.sh
```

## Troubleshooting

Memory Errors
```bash
# Error: STAR alignment fails due to insufficient memory
# Solution: Reduce STAR_JOBS or increase memory allocation
export STAR_JOBS=1  # Reduce parallel STAR jobs
# Or modify SLURM parameters: --mem-per-cpu=8G
```
Reference File Issues
```bash
# Error: Reference files not found
# Solution: Verify paths in conf/config.sh
ls -la /path/to/reference/files
# Update REFERENCE_DIR path accordingly
```

Parallel Processing Errors
```bash
# Error: Too many jobs running simultaneously
# Solution: Adjust parallelization parameters
export PARALLEL_JOBS=2    # Reduce from 8 to 2
export STAR_JOBS=1        # Reduce from 4 to 1
```

Log Files and Debugging

Important log files for troubleshooting:
```bash
# Main pipeline log
slurm/clinical_rna_Fusion_JOBID.out

# Step-specific logs
results/Input-trimmed/{sample}/{sample}.fastp.log
results/FastQC-trimmed/{sample}/{sample}.fastqc.log
results/Output/{sample}/Log.final.out  # STAR alignment log
results/Output/{sample}/arriba.log     # Arriba fusion log
```