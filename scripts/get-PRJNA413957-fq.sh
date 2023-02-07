#!/usr/bin/env bash
#
# Download HH1 10uM and DMSO 96hr fastq files from SRA
#
# ----------------------------------------------------------------------------
set -e
MEM=5G
THREADS=12

for RUN in $(tail -n+2 annotation-small.tsv | cut -f2) 
do
prefetch $RUN
fasterq-dump --mem $MEM --threads $THREADS $RUN
done

