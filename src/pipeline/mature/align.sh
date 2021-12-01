#!/bin/bash
set -e

[ -n "$SAMPLES" ] || { echo "SAMPLES is not set"; exit 1; }
[ -n "$BUILD_ROOT" ] || { echo "BUILD_ROOT is not set"; exit 1; }

mkdir -p "$BUILD_ROOT/mature/alignment"

for sample in $SAMPLES ; do
    filename="$BUILD_ROOT/short_long/${sample}_short_15-31.fq"
    bowtie -n 0 -l 32 --norc --best --strata -m 1 "$BUILD_ROOT/mature/index/mature" "$filename" \
     --un "$BUILD_ROOT/mature/alignment/${sample}_unaligned_maturemiRNA.fq" \
     -S "$BUILD_ROOT/mature/alignment/${sample}_aligned_maturemiRNA.sam" \
     2> "$BUILD_ROOT/mature/alignment/${sample}_maturemiRNA.log"
done
