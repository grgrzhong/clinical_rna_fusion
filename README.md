# Clinical RNA Fusion Analysis Pipeline

This pipeline performs RNA fusion detection using multiple tools including Arriba and STAR-Fusion. It processes paired-end RNA-seq FASTQ files through quality control, alignment, fusion detection, and generates comprehensive reports.

## Overview

The Clinical RNA Fusion Analysis Pipeline includes the following main steps:
1. **Preprocessing & QC**: FastQC analysis and read trimming with fastp
2. **Arriba Fusion Detection**: STAR alignment and Arriba fusion calling
3. **STAR-Fusion Detection**: STAR alignment and STAR-Fusion calling
4. **Fusion Report Generation**: Combine the Arriba and STAR fusion results
5. **QC Metrics**: Quality control metrics collection HS and RnaMetetics report

## Installation
### Step 1: Install Conda/Mamba

```bash
# Download and install Miniforge (includes conda and mamba)
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh

# Follow the installation prompts and restart your terminal
# Or source the conda initialization
source ~/.bashrc

# Verify installation
conda --version
```

### Step 2: Configure Conda Environment

```bash
# Update conda and configure channels
conda update -y -n base -c defaults conda
conda config --set solver libmamba
conda config --add channels conda-forge
conda config --set always_yes true
conda config --set auto_activate_base false
```

### Step 3: Create RNA Fusion Environment

```bash
# Create environment with required tools
conda create -n rnafusion apptainer parallel
```

### Step 4: Install System Dependencies (Optional)

```bash
# Install additional system dependencies for Apptainer (if needed)
sudo apt update
sudo apt install gocryptfs fuse2fs -y
```

### Step 5: Download the Pipeline

```bash
# Clone the repository
git clone https://github.com/grgrzhong/clinical_rna_fusion.git
cd clinical_rna_fusion

# Make the main pipeline script executable
chmod +x rnafusion_pipeline.sh
```

## Pipeline Steps

The pipeline consists of 5 main steps that are executed sequentially:

| Step | Script | Description |
|------|--------|-------------|
| 1 | `step_01_preprocess.sh` | Quality control (FastQC) and read trimming (fastp) |
| 2 | `step_02_arriba_fusion.sh` | STAR alignment and Arriba fusion detection |
| 3 | `step_03_star_fusion.sh` | STAR-Fusion based fusion detection |
| 4 | `step_04_fusion_report.sh` | Generate comprehensive fusion reports |
| 5 | `step_05_qc_metrics.sh` | Collect QC metrics and generate MultiQC report |

## Usage

### Basic Command Structure

```bash
./rnafusion_pipeline.sh <input_dir> <output_dir> <parallel_jobs> <star_jobs>
```

### Parameters

- **`input_dir`**: Absolute path to directory containing raw FASTQ files
- **`output_dir`**: Absolute path to directory for output files
- **`parallel_jobs`**: Number of samples to process in parallel
- **`star_jobs`**: Number of parallel STAR alignment jobs (each uses 16 threads)

> **Important**: Ensure `star_jobs × 16 ≤ total_CPU_cores` to avoid system overload.

### Example Commands

```bash
# Basic example with single sample processing
./rnafusion_pipeline.sh \
  /mnt/f/projects/clinical_rna_fusion/test/Raw \
  /mnt/f/projects/clinical_rna_fusion/test \
  1 1

# Processing multiple samples with more parallelization
./rnafusion_pipeline.sh \
  /path/to/fastq/files \
  /path/to/output \
  8 2
```

## Input Requirements

### FASTQ File Naming Convention

The pipeline expects paired-end FASTQ files with the following naming pattern:
```
{sample_id}_1.fastq.gz  # Forward reads
{sample_id}_2.fastq.gz  # Reverse reads
```

### Directory Structure Example

```
input_dir/
├── Sample1_1.fastq.gz
├── Sample1_2.fastq.gz
├── Sample2_1.fastq.gz
├── Sample2_2.fastq.gz
└── ...
```

## Output Structure

The pipeline generates the following output structure:

```
output_dir/
├── Input-trimmed/                              # Trimmed FASTQ files
├── FastQC-trimmed/                             # Trimmed FASTQ quality reports
├── Output/                                     # Main analysis results
    ├── {sample}/         
    │   ├── Aligned.sortedByCoord.out.bam       # Final bam file
    │   ├── fusions.tsv                         # Arriba fusion results
    │   ├── star-fusion.fusion_predictions.tsv  # STAR fusion results
    │   ├── {sample}_RnaSeqMetrics.csv          # Quality control metrics
    │   └── fusion_report/                      # combined Arriba and STAR fusion results
    
```