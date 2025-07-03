#!/bin/bash

INPUT=/media/hkusarc/Data/STUMP_RNA_seq_ALL/Input-trimmed   ##Modify if needed
REFERENCE=/media/hkusarc/Data/Reference/Gencode/gencode.hg38.v44/GRCh38.primary_assembly.genome.fa
ANNOTATION=/media/hkusarc/Data/Reference/Gencode/gencode.hg38.v44/gencode.v44.primary_assembly.annotation.gtf
BLACKLIST=/media/hkusarc/Data/Software/arriba_v2.4.0/database/blacklist_hg38_GRCh38_v2.4.0.tsv.gz
KNOWN_FUSION=/media/hkusarc/Data/Software/arriba_v2.4.0/database/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz
PROTEIN_DOMAINS=/media/hkusarc/Data/Software/arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.gff3
STARINDEX=/media/hkusarc/Data/Reference/Gencode/STAR_index_hg38.v44
CYTOBANDS=/media/hkusarc/Data/Software/arriba_v2.4.0/database/cytobands_hg38_GRCh38_v2.4.0.tsv

export PATH="/media/hkusarc/Data/Software/arriba_v2.4.0:/media/hkusarc/Data/Software/STAR-2.7.11b/bin/Linux_x86_64:$PATH"



for file in $(ls ${INPUT}/*.fastq.gz | grep "R1"); do 
echo $file; 
FILENAME=$(basename "$file" | cut -d "_" -f 1)
echo "Processing file: $file"
echo "Sample name: $FILENAME"

OUTPUT="${INPUT}/../Output_Arriba/${FILENAME}/"
mkdir -p "$OUTPUT"  ## Create the directory if it doesn't exist

# STAR \
#     --runThreadN 16 \
#     --genomeDir $STARINDEX \
#     --readFilesIn ${file} ${file//R1/R2} \
#     --readFilesCommand "gunzip -c" \
#     --outFileNamePrefix $OUTPUT \
#     --outSAMtype BAM SortedByCoordinate \
#     --outSAMunmapped Within \
#     --outBAMcompression 0 \
#     --outFilterMultimapNmax 50 \
#     --peOverlapNbasesMin 10 \
#     --alignSplicedMateMapLminOverLmate 0.5 \
#     --alignSJstitchMismatchNmax 5 -1 5 5 \
#     --chimSegmentMin 10 \
#     --chimOutType WithinBAM HardClip \
#     --chimJunctionOverhangMin 10 \
#     --chimScoreDropMax 30 \
#     --chimScoreJunctionNonGTAG 0 \
#     --chimScoreSeparation 1 \
#     --chimSegmentReadGapMax 3 \
#     --chimMultimapNmax 50;
   
# #    --twopassMode Basic \

# samtools index ${OUTPUT}/Aligned.sortedByCoord.out.bam
  
arriba \
    -x ${OUTPUT}/Aligned.sortedByCoord.out.bam \
    -o ${OUTPUT}/fusions.tsv \
    -O ${OUTPUT}/fusions.discarded.tsv \
    -a $REFERENCE -g $ANNOTATION \
    -b $BLACKLIST -k $KNOWN_FUSION -t $KNOWN_FUSION -p $PROTEIN_DOMAINS;


/media/hkusarc/Data/Software/arriba_v2.4.0/draw_fusions.R \
	--fusions=${OUTPUT}/fusions.tsv \
	--output=${OUTPUT}/fusions.pdf \
	--annotation=$ANNOTATION \
	--cytobands=$CYTOBANDS \
	--proteinDomains=$PROTEIN_DOMAINS \
	--minConfidenceForCircosPlot=none \
	--mergeDomainsOverlappingBy=0.5
       

done

##--chimSegmentMin default: 0; int>=0: minimum length of chimeric segment length, if ==0, no chimeric output

##--chimJunctionOverhangMin default: 20; int>=0: minimum overhang for a chimeric junction

##--chimOutJunctionFormat; default: 0; int: formatting type for the Chimeric.out.junction file; 
##0 no comment lines/headers 1 comment lines at the end of the file: command line and Nreads: total, unique/multi-mapping

##--alignSJDBoverhangMin default: 3 int>0: minimum overhang (i.e. block size) for annotated (sjdb) spliced alignments

##--alignSJoverhangMin default: 5 int>0: minimum overhang (i.e. block size) for spliced alignments

##--alignMatesGapMax default: 0 maximum gap between two mates, if 0, max intron gap will be determined by (2ˆwinBinNbits)*winAnchorDistNbins

##--alignIntronMax default: 0 maximum intron size, if 0, max intron size will be determined by (2ˆwinBinNbits)*winAnchorDistNbins

##--alignSJstitchMismatchNmax default: 0 -1 0 0; 4*int>=0: maximum number of mismatches for stitching of the splice junctions (-1: no limit). (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif.

## --chimMultimapScoreRange default: 1 int>=0: the score range for multi-mapping chimeras below the best chimeric score. Only works with –chimMultimapNmax > 1

## --chimScoreJunctionNonGTAG default: -1 int: penalty for a non-GT/AG chimeric junction

##--chimMultimapNmax default: 0 int>=0: maximum number of chimeric multi-alignments 0 use the old scheme for chimeric detection which only considered unique alignments

## --chimNonchimScoreDropMin default: 20 52int>=0: to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be greater than this value


### --outSAMstrandField intronMotif \ # include for potential use with StringTie for assembly
### --chimSegmentMin 12 \ # ** essential to invoke chimeric read detection & reporting **
### --chimOutJunctionFormat 1  \ # **essential** includes required metadata in Chimeric.junction.out file.
### --alignMatesGapMax 100000 \ # avoid readthru fusions within 100k
### --alignSJstitchMismatchNmax 5 -1 5 5 \  # settings improved certain chimera detections
