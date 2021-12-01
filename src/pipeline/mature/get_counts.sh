#!/bin/bash
set -e

[ -n "$SAMPLES" ] || { echo "SAMPLES is not set"; exit 1; }
[ -n "$BUILD_ROOT" ] || { echo "BUILD_ROOT is not set"; exit 1; }

mkdir -p "$BUILD_ROOT/mature/counts"

for sample in $SAMPLES ; do
    samtools idxstats "$BUILD_ROOT/mature/alignment/${sample}_aligned_maturemiRNA.bam" | cut -f1,3 - | \
     sed -e "1s/^/miRNA\t${sample}-miRNAcount\n/" > "$BUILD_ROOT/mature/counts/${sample}_mature_counts.tsv"
done
