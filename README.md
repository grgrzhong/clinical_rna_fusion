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

## Create environment with required tools
conda create -n rnafusion apptainer parallel

## Install apptainer required software
sudo apt install gocryptfs -y
sudo apt install fuse2fs -y

```

## Runing pipeline

Arguments to run the rnafusion pipeline:

- input_dir: directory contains raw fastq files, must be absolute path
- output_dir: directory to save the output files, must be absolute path
- parallel_jobs: number of jobs in parallel procesing samples
- star_jobs: number of jobs in parallel parallel procesing samples for running STAR alignment, each job takes 16 threads (default)

```bash
rnafusion_pipeline <input_dir> <output_dir> <parallel_jobs> <star_jobs>
```