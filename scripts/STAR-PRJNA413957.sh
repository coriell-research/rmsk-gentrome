#!/usr/bin/env bash
#
# Align reads with STAR for input into TEtranscripts
#
# ----------------------------------------------------------------------------
FQ=/home/gcalendo/data/rmsk-gentrome/data/PRJNA413957
SAMPLES=/home/gcalendo/data/rmsk-gentrome/scripts/sample-names.tsv
IDX=/mnt/data/gdata/human/hg38/GENCODE/STAR_idx/
OUT=/home/gcalendo/data/rmsk-gentrome/data/PRJNA413957/STAR_outs
MM=100
THREADS=32

while read -r SAMPLE FQ1 FQ2
do
STAR --runThreadN $THREADS \
     --genomeDir $IDX \
     --readFilesIn ${FQ}/${FQ1},${FQ}/${FQ2} \
     --readFilesCommand zcat \
     --outFilterType BySJout \
     --outFileNamePrefix ${OUT}/${SAMPLE}. \
     --outFilterMultimapNmax $MM \
     --winAnchorMultimapNmax $MM \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverReadLmax 0.04 \
     --alignIntronMin 20 \
     --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000 \
     --outMultimapperOrder Random \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts;
done < $SAMPLES

