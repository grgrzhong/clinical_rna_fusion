#!/bin/bash

## Generate and Download containers from the Seqera community registry using Singularity/Apptainer. 
## https://seqera.io/containers/

## Make sure you have singularity or apptainer installed on your system.
container_dir="${PWD}/containers"
mkdir -p "${container_dir}"

singularity pull --force --dir "${container_dir}" arriba.sif oras://community.wave.seqera.io/library/arriba:2.4.0--d895079af70ac3a9
singularity pull --force --dir "${container_dir}" fusion_report.sif oras://community.wave.seqera.io/library/fusion-report:4.0.1--0e79b958dd0158e0
singularity pull --force --dir "${container_dir}" picard.sif oras://community.wave.seqera.io/library/picard:3.4.0--2976616e7cbd4840

## Download STAR-Fusion containers https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/STAR-Fusion/
singularity pull --force --dir ${container_dir} star-fusion.v1.15.0.simg https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/STAR-Fusion/star-fusion.v1.15.0.simg
singularity pull --force --dir ${container_dir} star_fusion.sif https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/STAR-Fusion/star-fusion.v1.15.1.simg