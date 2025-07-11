#!/bin/bash

# Move the outputs to multiomics server

multiomics_dir=/mnt/m/RNA-seq/Exome

## Copy the outputs in the output foder to the output folder in multiomics server
rsync -av --progress  --no-g --no-t \
    $PWD/data/Output/ \
    $multiomics_dir/Output/

rsync -av --progress  --no-g --no-t \
    $PWD/data/Input-trimmed/ \
    $multiomics_dir/Input-trimmed/

rsync -av --progress  --no-g --no-t \
    $PWD/data/FastQC-trimmed/ \
    $multiomics_dir/FastQC-trimmed/
