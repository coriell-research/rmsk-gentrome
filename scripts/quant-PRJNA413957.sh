#!/usr/bin/env bash
#
# Quantify reads with Salmon using the new transcriptome index
#
# NOTE: The samples are split across Runs, they are NOT paired-end
# -----------------------------------------------------------------------------
set -e
SAMPLES=/home/gcalendo/data/rmsk-gentrome/scripts/sample-names.tsv
IDX=/home/gcalendo/data/rmsk-gentrome/data/GRCh38-rmsk-salmon-idx
FQ=/home/gcalendo/data/rmsk-gentrome/data/PRJNA413957
OUT=${FQ}/quants
THREADS=32
BOOT=30

while read -r SAMPLE FQ1 FQ2
do
salmon quant \
  -i $IDX \
  -l A \
  -r $FQ/$FQ1 $FQ/$FQ2 \
  --validateMappings \
  --gcBias \
  --seqBias \
  --threads $THREADS \
  --numBootstraps $BOOT \
  -o $OUT/${SAMPLE}_quants
done < $SAMPLES

