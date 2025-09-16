#!/bin/bash

# sbatch /lustre1/g/path_my/pipeline/clinical_rna_fusion/rnafusion_pipeline.sh \
#     /lustre1/g/path_my/pipeline/clinical_rna_fusion/data/DFSP/YeungMCF_UCSO_CPOS-250731-HTS-30262a/primary_seq \
#     /lustre1/g/path_my/pipeline/clinical_rna_fusion/results/DFSP/30262a \
#     16 \
#     2 \
#     /lustre1/g/path_my/Reference

# sbatch /lustre1/g/path_my/pipeline/clinical_rna_fusion/rnafusion_pipeline.sh \
#     /lustre1/g/path_my/pipeline/clinical_rna_fusion/data/DFSP/YeungMCF_UCSO_CPOS-250808-HTS-30362a/primary_seq \
#     /lustre1/g/path_my/pipeline/clinical_rna_fusion/results/DFSP/30362a \
#     16 \
#     2 \
#     /lustre1/g/path_my/Reference

sbatch /lustre1/g/path_my/pipeline/clinical_rna_fusion/rnafusion_pipeline.sh \
    /lustre1/g/path_my/pipeline/clinical_rna_fusion/data/DFSP/YeungMCF_UCSO_CPOS-250815-HTS-30439a/primary_seq \
    /lustre1/g/path_my/pipeline/clinical_rna_fusion/results/DFSP/30439a \
    16 \
    2 \
    /lustre1/g/path_my/Reference