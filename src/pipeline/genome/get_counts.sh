#!/bin/bash
set -e

[ -n "$SAMPLES" ] || { echo "SAMPLES is not set"; exit 1; }
[ -n "$BUILD_ROOT" ] || { echo "BUILD_ROOT is not set"; exit 1; }
[ -n "$DATA_ROOT" ] || { echo "$DATA_ROOT is not set"; exit 1; }
[ -n "$SRC_ROOT" ] || { echo "SRC_ROOT is not set"; exit 1; }

mkdir -p "$BUILD_ROOT/genome/counts"

for sample in $SAMPLES ; do
  filename="$BUILD_ROOT/genome/alignment/${sample}_aligned_genome_tagged.bam"
  printf "\nProcessing sample: %s\n" "$sample"
  samtools view "$filename" | \
    "$SRC_ROOT/pipeline/utils/get_counts.py" -m "$DATA_ROOT/mirnas.txt" -o "$BUILD_ROOT/genome/counts/${sample}_genome_based_counts.txt"
  sed "N;s/\n/\t/g" "$BUILD_ROOT/genome/counts/${sample}_genome_based_counts.txt" | \
   sort -k1 | uniq | sed "1s/^/miRNA\t${sample}\n/" > "$BUILD_ROOT/genome/counts/${sample}_taggedBAMcounts.txt"
done