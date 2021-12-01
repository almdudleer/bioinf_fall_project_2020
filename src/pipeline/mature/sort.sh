#!/bin/bash
set -e

[ -n "$SAMPLES" ] || { echo "SAMPLES is not set"; exit 1; }
[ -n "$BUILD_ROOT" ] || { echo "BUILD_ROOT is not set"; exit 1; }

for sample in $SAMPLES ; do
    samtools sort "$BUILD_ROOT/mature/alignment/${sample}_aligned_maturemiRNA.sam" > "$BUILD_ROOT/mature/alignment/${sample}_aligned_maturemiRNA.bam"
    samtools index "$BUILD_ROOT/mature/alignment/${sample}_aligned_maturemiRNA.bam"
done
