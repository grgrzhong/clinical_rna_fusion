#!/bin/bash

## Generate and Download containers from the Seqera community registry using Singularity/Apptainer. 
## https://seqera.io/containers/

## Make sure you have singularity or apptainer installed on your system.

singularity pull --force --dir "${CONTAINER_DIR}" arriba.sif oras://community.wave.seqera.io/library/arriba:2.4.0--d895079af70ac3a9
singularity pull --force --dir "${CONTAINER_DIR}" fastp.sif oras://community.wave.seqera.io/library/fastp:0.24.0--0397de619771c7ae
singularity pull --force --dir "${CONTAINER_DIR}" fastqc.sif oras://community.wave.seqera.io/library/fastqc:0.12.1--104d26ddd9519960
singularity pull --force --dir "${CONTAINER_DIR}" fusion_report.sif oras://community.wave.seqera.io/library/fusion-report:4.0.1--0e79b958dd0158e0
singularity pull --force --dir "${CONTAINER_DIR}" multiqc.sif oras://community.wave.seqera.io/library/multiqc:1.28--d466e41d58d6d704
singularity pull --force --dir "${CONTAINER_DIR}" picard.sif oras://community.wave.seqera.io/library/picard:3.4.0--2976616e7cbd4840
singularity pull --force --dir "${CONTAINER_DIR}" r.sif oras://community.wave.seqera.io/library/r-base_r-fs_r-here_r-optparse_pruned:6d3c4357c207ae65
singularity pull --force --dir "${CONTAINER_DIR}" samtools.sif oras://community.wave.seqera.io/library/samtools:1.21--84c9d77c3901e90b

## Download STAR-Fusion containers https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/STAR-Fusion/
singularity pull --force --dir "${CONTAINER_DIR}" star-fusion.sif https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/STAR-Fusion/star-fusion.v1.15.1.simg

singularity pull --force --dir "${CONTAINER_DIR}" subread.sif oras://community.wave.seqera.io/library/subread:2.1.1--bae420bffb4edf16
singularity pull --force --dir "${CONTAINER_DIR}" rsem.sif library://yh549848/rnaseqde/rsem:1.3.3
## Download RSEM container for transcript quantification
# singularity pull --force --dir "${CONTAINER_DIR}" rsem.sif oras://community.wave.seqera.io/library/rsem:1.3.3--3b6cfedbe5b2229e
# singularity pull --arch amd64 library://yh549848/rnaseqde/rsem:1.3.3
# singularity pull --arch amd64 library://yh549848/rnaseqde/rsem:sha256.d2f1b7a0c1c94db46cf6d3ca53a620975c2932fd99371198de1db509541701bd