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

## Create environment with required software
conda create -n rnafusion 

```

## Runing pipeline
```bash
rnafusion_pipeline <primary_seq_dir> <output_dir> <parallel_jobs>

## example
rnafusion_pipeline 
```