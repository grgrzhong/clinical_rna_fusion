## Install software
```bash
## Install conda
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-$(uname)-$(uname -m).sh

## Add conda to PATH if found
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    . $HOME/miniconda3/etc/profile.d/conda.sh
fi

## General config
conda update -y -n base -c defaults conda
conda config --set solver libmamba
conda config --add channels conda-forge
conda config --set always_yes true
conda config --set auto_activate_base false

## Create environment with required software
conda create -n rnafusion apptainer parallel

```

## Runing pipeline

- primary_seq_dir: directory contains raw fastq files
- output_dir: directory to save the output files
- parallel_jobs: number of jobs in parallel procesing samples

```bash
rnafusion_pipeline <primary_seq_dir> <output_dir> <parallel_jobs>

## example
rnafusion_pipeline 
```