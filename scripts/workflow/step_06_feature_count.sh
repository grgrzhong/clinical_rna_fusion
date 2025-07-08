#!/bin/bash
ANNOTATION=/media/maximus/Data/Reference/Gencode/gencode.v36.primary_assembly.annotation.gtf

samtools sort -o /media/maximus/Data/RNA-seq/Output/Test/Aligned.sortedByCoord.out.bam

featureCounts  \
		-T 4 \
		-s 2 \
		-a $ANNOTATION \
		-o /media/maximus/Data/RNA-seq/Output/Test/Test_featureCounts.txt \
		/media/maximus/Data/RNA-seq/Output/Test/Aligned.sortedByCoord.out.bam

		
