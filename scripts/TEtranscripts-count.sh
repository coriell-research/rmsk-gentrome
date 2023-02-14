#!/usr/bin/env bash
#
# Count aligned reads to TEs with TEtranscripts
#
# ----------------------------------------------------------------------------
GTF=/mnt/data/gdata/human/hg38/GENCODE/TEtranscripts/gencode.v38.annotation.gtf.ind
TE=/mnt/data/gdata/human/hg38/GENCODE/TEtranscripts/GRCh38_GENCODE_rmsk_TE.gtf.ind
BAM=/home/gcalendo/data/rmsk-gentrome/data/PRJNA413957/STAR_outs
SAMPLES=biosample-names.txt
OUT=/home/gcalendo/data/rmsk-gentrome/data/PRJNA413957/TEtranscripts_outs
JOBS=6
LEN=50

parallel --jobs $JOBS \
  "TEcount -b $BAM/{}.Aligned.sortedByCoord.out.bam \
           --GTF $GTF \
           --TE $TE \
           --format BAM \
           --mode multi \
           --project {} \
           --outdir $OUT \
           --stranded reverse \
           --sortByPos \
           --fragmentLength $LEN" :::: $SAMPLES
 
