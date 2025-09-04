# Clinical RNA Fusion Analysis Pipeline

![Pipeline Status](https://img.shields.io/badge/pipeline-active-brightgreen.svg)
![Platform](https://img.shields.io/badge/platform-Linux-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

A simple pipeline for detecting RNA fusions in clinical samples using multiple state-of-the-art tools including Arriba and STAR-Fusion. This pipeline processes paired-end RNA-seq FASTQ files through quality control, alignment, fusion detection, and generates detailed reports suitable for clinical applications.

## 🔬 Overview

The Clinical RNA Fusion Analysis Pipeline is designed for robust and reproducible fusion detection in clinical RNA-seq data. The workflow includes:

1. **Preprocessing & QC**: Quality assessment with FastQC and adapter trimming with fastp
2. **Arriba Fusion Detection**: STAR alignment followed by Arriba fusion calling with blacklist filtering
3. **STAR-Fusion Detection**: Dedicated STAR-Fusion analysis with CTAT genome library
4. **Fusion Report Generation**: Comprehensive fusion report combining results from multiple tools
5. **Quality Control Metrics**: RNA-seq and hybrid selection metrics collection with MultiQC reporting

## 🚀 Features

- **Multi-tool fusion detection**: Combines Arriba and STAR-Fusion for enhanced sensitivity and specificity
- **Containerized execution**: Uses Singularity/Apptainer containers for reproducibility
- **Parallel processing**: Configurable parallelization for efficient resource utilization
- **Comprehensive QC**: Multiple quality control checkpoints and metrics
- **Clinical-ready**: Includes known fusion databases and blacklist filtering
- **HPC compatible**: SLURM batch system integration for high-performance computing
- **Flexible configuration**: Command-line parameters and configuration files

## 📋 Prerequisites

### System Requirements
- **OS**: Linux (Ubuntu 18.04+ or CentOS 7+ recommended)
- **Memory**: Minimum 32GB RAM (64GB+ recommended for multiple samples)
- **Storage**: 100GB+ free space for references and intermediate files
- **CPU**: Multi-core processor (16+ cores recommended)

### Software Dependencies
- **Conda/Mamba**: Package management
- **Singularity/Apptainer**: Container runtime
- **GNU Parallel**: Parallel job execution

## 🛠 Installation

### Step 1: Install Conda/Mamba

```bash
# Download and install Miniforge (includes conda and mamba)
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh

# Follow installation prompts and restart terminal
source ~/.bashrc

# Verify installation
conda --version
mamba --version
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
# Create environment from configuration file
conda env create -f ./conf/rnafusion.yml

# Or create manually with required tools
conda create -n rnafusion apptainer parallel r-base -c conda-forge
```

### Step 4: Setup Pipeline

```bash
# Make scripts executable
chmod +x rnafusion_pipeline.sh
chmod +x run_pipeline_local.sh
chmod +x run_pipeline_hpc.sh
chmod +x scripts/workflow/*.sh
```

## 🔧 Configuration

### Environment Variables

The pipeline uses `conf/config.sh` for configuration. Key parameters:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `INPUT_DIR` | Raw FASTQ files directory | `${PROJECT_DIR}/data/test_data/Raw` |
| `OUTPUT_DIR` | Results output directory | `${PROJECT_DIR}/results` |
| `PARALLEL_JOBS` | Parallel sample processing | `1` |
| `STAR_JOBS` | STAR alignment jobs (16 threads each) | `1` |
| `REFERENCE_DIR` | Reference files directory | `/mnt/f/Reference` |
| `CONTAINER_DIR` | Singularity containers directory | `${PROJECT_DIR}/containers` |

### Resource Allocation

> ⚠️ **Important**: Ensure `STAR_JOBS × 16 ≤ total_CPU_cores` to prevent system overload.

For optimal performance:
- **Single sample**: `PARALLEL_JOBS=1, STAR_JOBS=1` (16 cores)
- **Multiple samples (32 cores)**: `PARALLEL_JOBS=4, STAR_JOBS=2` 
- **High-throughput (64+ cores)**: `PARALLEL_JOBS=8, STAR_JOBS=4`

## 📊 Usage

### Local Execution

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

### HPC/SLURM Execution

```bash
# Submit to SLURM scheduler
./run_pipeline_hpc.sh

# Or submit directly
sbatch rnafusion_pipeline.sh
```

### Parameter Options

```bash
./rnafusion_pipeline.sh [INPUT_DIR] [OUTPUT_DIR] [PARALLEL_JOBS] [STAR_JOBS] [REFERENCE_DIR] [CONTAINER_DIR]
```

**Examples:**

```bash
# Single sample processing
./rnafusion_pipeline.sh \
    /data/samples/batch1 \
    /results/batch1 \
    1 1

# Multi-sample processing
./rnafusion_pipeline.sh \
    /data/samples/batch2 \
    /results/batch2 \
    4 \
    2 \
    /references/hg38 \
    /containers
```

## 📁 Input Requirements

### FASTQ File Naming Convention

The pipeline expects paired-end FASTQ files following this pattern:
```
{sample_id}_1.fastq.gz  # Forward reads (R1)
{sample_id}_2.fastq.gz  # Reverse reads (R2)
```

### Example Input Structure

```
input_directory/
├── Patient001_1.fastq.gz
├── Patient001_2.fastq.gz
├── Patient002_1.fastq.gz
├── Patient002_2.fastq.gz
├── Control001_1.fastq.gz
├── Control001_2.fastq.gz
└── ...
```

### File Requirements
- **Format**: Gzipped FASTQ (.fastq.gz)
- **Read type**: Paired-end sequencing
- **Quality**: Phred33 encoding
- **Size**: No specific limits (tested with 1GB-50GB files)

## 📈 Output Structure

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

## 🔍 Pipeline Steps

| Step | Script | Tools Used | Description |
|------|--------|------------|-------------|
| **1** | `step_01_preprocess.sh` | fastp, FastQC | Quality control and adapter trimming |
| **2** | `step_02_arriba_fusion.sh` | STAR, Arriba | STAR alignment and Arriba fusion detection |
| **3** | `step_03_star_fusion.sh` | STAR-Fusion | Dedicated STAR-Fusion analysis |
| **4** | `step_04_fusion_report.sh` | Custom R scripts | Generate comprehensive fusion reports |
| **5** | `step_05_qc_metrics.sh` | Picard, MultiQC | Collect QC metrics and generate reports |

### Detailed Step Information

#### Step 1: Preprocessing
- **fastp**: Adapter trimming, quality filtering, poly-G/X trimming
- **FastQC**: Quality assessment of trimmed reads
- **Output**: Clean reads ready for alignment

#### Step 2: Arriba Fusion Detection
- **STAR**: Splice-aware alignment with chimeric read detection
- **Arriba**: Fusion detection with blacklist filtering
- **Features**: Known fusion annotation, protein domain analysis

#### Step 3: STAR-Fusion Detection
- **STAR-Fusion**: Alternative fusion detection method
- **CTAT Library**: Comprehensive cancer transcriptome resource
- **Filtering**: Multiple evidence requirements and confidence scoring

#### Step 4: Fusion Reporting
- **Integration**: Combines results from Arriba and STAR-Fusion
- **Annotation**: Adds gene information and clinical relevance
- **Formats**: Excel and JSON outputs for downstream analysis

#### Step 5: Quality Control
- **Picard**: RNA-seq and hybrid selection metrics
- **MultiQC**: Integrated quality report across all samples
- **Metrics**: Read alignment, insert size, coverage statistics

## 🧬 Containers and Dependencies

The pipeline uses Singularity/Apptainer containers for reproducibility:

| Tool | Container | Version | Purpose |
|------|-----------|---------|---------|
| fastp | `fastp.sif` | Latest | Read trimming and QC |
| FastQC | `fastqc.sif` | Latest | Quality assessment |
| STAR | `star-fusion.sif` | 1.15.0 | Read alignment |
| Arriba | `arriba.sif` | 2.4.0 | Fusion detection |
| STAR-Fusion | `star-fusion.sif` | 1.15.0 | Alternative fusion detection |
| Picard | `picard.sif` | Latest | QC metrics collection |
| MultiQC | `multiqc.sif` | Latest | Report aggregation |
| R | `r.sif` | Latest | Statistical analysis |
| Samtools | `samtools.sif` | Latest | BAM file processing |

Containers are automatically downloaded and configured by the `conf/containers.sh` script.

## 📋 Quality Control Metrics

The pipeline generates comprehensive quality metrics:

### RNA-seq Metrics (Picard CollectRnaSeqMetrics)
- **Alignment metrics**: Total reads, aligned reads, mapping rate
- **RNA-specific metrics**: Ribosomal RNA fraction, coding/UTR/intronic rates
- **Insert size distribution**: Fragment length statistics
- **Strand specificity**: Sense/antisense mapping ratios

### Hybrid Selection Metrics (Picard CollectHsMetrics)
- **Target coverage**: On-target vs off-target reads
- **Coverage uniformity**: Target coverage distribution
- **Enrichment efficiency**: Fold enrichment over genome
- **Coverage depth**: Mean/median target coverage

### Fusion-specific Metrics
- **Detection sensitivity**: Number of fusion candidates
- **Tool concordance**: Overlap between Arriba and STAR-Fusion
- **Confidence scoring**: High/medium/low confidence fusions
- **Known fusion annotation**: Clinically relevant fusion matches

## 🚨 Troubleshooting

### Common Issues

#### Memory Errors
```bash
# Error: STAR alignment fails due to insufficient memory
# Solution: Reduce STAR_JOBS or increase memory allocation
export STAR_JOBS=1  # Reduce parallel STAR jobs
# Or modify SLURM parameters: --mem-per-cpu=8G
```
#### Reference File Issues
```bash
# Error: Reference files not found
# Solution: Verify paths in conf/config.sh
ls -la /path/to/reference/files
# Update REFERENCE_DIR path accordingly
```

#### Parallel Processing Errors
```bash
# Error: Too many jobs running simultaneously
# Solution: Adjust parallelization parameters
export PARALLEL_JOBS=2    # Reduce from 8 to 2
export STAR_JOBS=1        # Reduce from 4 to 1
```

### Log Files and Debugging

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

### Performance Optimization

#### For Large Cohorts (50+ samples)
```bash
# Use higher parallelization
./rnafusion_pipeline.sh INPUT OUTPUT 8 4

# Consider batch processing
# Process samples in groups of 20-30
```

#### For Limited Resources
```bash
# Sequential processing
./rnafusion_pipeline.sh INPUT OUTPUT 1 1

# Use smaller reference indices if available
```