# Clinical RNA Fusion Pipeline

A pipeline for detecting RNA fusions in clinical samples using Arriba and STAR-Fusion.

## Quick Start

```bash
# Run with default settings
./run_pipeline_hpc.sh

# Run with custom parameters
./rnafusion_pipeline.sh <input_dir> <output_dir> <parallel_jobs> <star_jobs> <reference_dir> <container_dir>
```

## Setup

```bash
# Install Miniforge
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh
source ~/.bashrc

# Create environment
conda env create -f ./conf/base.yml
```

## Input
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
## Output

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
./run_pipeline_hpc.sh

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