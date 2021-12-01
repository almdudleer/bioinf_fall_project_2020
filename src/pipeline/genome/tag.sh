#!/bin/bash
set -e

[ -n "$SAMPLES" ] || { echo "SAMPLES is not set"; exit 1; }
[ -n "$BUILD_ROOT" ] || { echo "BUILD_ROOT is not set"; exit 1; }
[ -n "$DATA_ROOT" ] || { echo "DATA_ROOT is not set"; exit 1; }

for sample in $SAMPLES ; do
  bedtools tag -i "$BUILD_ROOT/genome/alignment/${sample}_aligned_genome_renamed.bam" \
    -files "$DATA_ROOT/hsa-genome-miRBase22v-onlymiRNAs-convforTagBAM.bed" -names -tag XQ \
    > "$BUILD_ROOT/genome/alignment/${sample}_aligned_genome_tagged.bam"
done
