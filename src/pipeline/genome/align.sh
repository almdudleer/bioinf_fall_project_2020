#!/bin/bash
set -e

[ -n "$SAMPLES" ] || { echo "SAMPLES is not set"; exit 1; }
[ -n "$BUILD_ROOT" ] || { echo "BUILD_ROOT is not set"; exit 1; }

mkdir -p "$BUILD_ROOT/genome/alignment"

for sample in $SAMPLES ; do
  filename="$BUILD_ROOT/mature/alignment/${sample}_unaligned_maturemiRNA.fq"
  bowtie -n 1 -l 32 --norc --best --strata -m 1 \
   "$BUILD_ROOT/genome/index/GCA_000001405.15_GRCh38_no_alt_analysis_set" "$filename" \
   --al "$BUILD_ROOT/genome/alignment/${sample}_aligned_genome.fastq" \
   -S "$BUILD_ROOT/genome/alignment/${sample}_aligned_genome.sam" \
   2> "$BUILD_ROOT/genome/alignment/${sample}_aligned_genomeaftermiRNA.log"
done
