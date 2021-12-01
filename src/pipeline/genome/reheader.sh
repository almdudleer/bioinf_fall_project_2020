#!/bin/bash
set -e

[ -n "$SAMPLES" ] || { echo "SAMPLES is not set"; exit 1; }
[ -n "$BUILD_ROOT" ] || { echo "BUILD_ROOT is not set"; exit 1; }
[ -n "$DATA_ROOT" ] || { echo "DATA_ROOT is not set"; exit 1; }

for sample in $SAMPLES ; do
    header=$(samtools view -H "$BUILD_ROOT/genome/alignment/${sample}_aligned_genome.bam")
    map="$DATA_ROOT/ncbitochrom.txt"
    while IFS=$'\t' read -r ucsc_name ncbi_name; do
        header=${header//$ncbi_name/$ucsc_name/}
    done < "$map"
    echo "$header" > "$BUILD_ROOT/genome/alignment/${sample}_tmp_header.txt"
    samtools reheader "$BUILD_ROOT/genome/alignment/${sample}_tmp_header.txt" \
      "$BUILD_ROOT/genome/alignment/${sample}_aligned_genome.bam" \
      > "$BUILD_ROOT/genome/alignment/${sample}_aligned_genome_renamed.bam"
    rm "$BUILD_ROOT/genome/alignment/${sample}_tmp_header.txt"
done
