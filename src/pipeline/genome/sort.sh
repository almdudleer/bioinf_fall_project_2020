#!/bin/bash
set -e

[ -n "$SAMPLES" ] || { echo "SAMPLES is not set"; exit 1; }
[ -n "$BUILD_ROOT" ] || { echo "BUILD_ROOT is not set"; exit 1; }

for sample in $SAMPLES ; do
  filename="$BUILD_ROOT/genome/alignment/${sample}_aligned_genome.sam"
  samtools sort "$filename" > "$BUILD_ROOT/genome/alignment/${sample}_aligned_genome.bam"
  samtools index "$BUILD_ROOT/genome/alignment/${sample}_aligned_genome.bam"
done
