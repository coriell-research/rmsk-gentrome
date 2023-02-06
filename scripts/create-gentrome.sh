#!/usr/bin/env bash
#
# Create the Salmon index for the rmsk-gentrome
#
# ----------------------------------------------------------------------------
FA=/home/gcalendo/data/rmsk-gentrome/data/rmsk-gentrome.fa.gz
DECOYS=/home/gcalendo/data/rmsk-gentrome/data/decoys.txt
IDX=/home/gcalendo/data/rmsk-gentrome/data/GRCh38-rmsk-salmon-idx
THREADS=24

salmon index --transcripts $FA -i $IDX --threads $THREADS --decoys $DECOYS --gencode

