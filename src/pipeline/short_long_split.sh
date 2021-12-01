#!/bin/bash
set -e

[ -n "$SAMPLES" ] || { echo "SAMPLES is not set"; exit 1; }
[ -n "$BUILD_ROOT" ] || { echo "BUILD_ROOT is not set"; exit 1; }

mkdir -p "$BUILD_ROOT/short_long"

for sample in $SAMPLES ; do
  filename="$BUILD_ROOT/revcomp/${sample}_trim_revcomp.fq.gz"
  seqkit seq "$filename" -M 31 > "$BUILD_ROOT/short_long/${sample}_short.fq"
  seqkit seq "$filename" -m 32 > "$BUILD_ROOT/short_long/${sample}_long.fq"
  seqkit seq "$BUILD_ROOT/short_long/${sample}_short.fq" -m 15 > "$BUILD_ROOT/short_long/${sample}_short_15-31.fq"
done
