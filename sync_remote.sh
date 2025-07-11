#!/bin/bash

# Move the outputs to multiomics server

local_path=$(pwd)
remote_path=/mnt/m/RNA-seq/Exome

## Copy the outputs in the output foder to the output folder in multiomics server
rsync -av --update --progress  --no-g --no-t \
    $local_path/data/Output/ \
    $remote_path/Output/

rsync -av --update --progress  --no-g --no-t \
    $local_path/data/Input-trimmed/ \
    $remote_path/Input-trimmed/

rsync -av --update --progress  --no-g --no-t \
    $local_path/data/FastQC-trimmed/ \
    $remote_path/FastQC-trimmed/
