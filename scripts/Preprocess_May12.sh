#!/bin/bash

INPUT=/media/hkusarc/Data/clinical_rnaseq_jun30   ##Modify if needed
OUTPUT=$(realpath "${INPUT}/../Input-trimmed") 

if [ ! -d "$OUTPUT" ]; then
  mkdir -p "$OUTPUT"
fi

for file in $(ls ${INPUT}/*.fastq.gz | rev | cut -d '_' -f 2- | rev | uniq); do
FILENAME=`echo $(basename $file)`; 
echo $FILENAME;

/media/hkusarc/Data/Software/fastp \
-i ${file}_1.fastq.gz \
-I ${file}_2.fastq.gz \
-o ${OUTPUT}/${FILENAME}_trimmed_R1.fastq.gz \
-O ${OUTPUT}/${FILENAME}_trimmed_R2.fastq.gz \
--detect_adapter_for_pe \
--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
--adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
--dont_eval_duplication \
--trim_poly_g \
--trim_poly_x \
-j ${OUTPUT}/${FILENAME}.json \
-h ${OUTPUT}/${FILENAME}.html \
-w 8

done

mkdir ${OUTPUT}/FastQC_trimmed

/media/hkusarc/Data/Software/FastQC/fastqc \
${OUTPUT}/*.fastq.gz \
-o ${OUTPUT}/FastQC_trimmed \
--memory 8192 -t 8 


