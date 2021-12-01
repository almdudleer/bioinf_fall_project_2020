#!/bin/bash
set -e

[ -n "$SAMPLES" ] || { echo "SAMPLES is not set"; exit 1; }
[ -n "$BUILD_ROOT" ] || { echo "BUILD_ROOT is not set"; exit 1; }

mkdir -p "$BUILD_ROOT/revcomp"

for sample in $SAMPLES ; do
  seqkit seq -r -p "$BUILD_ROOT/trimmed/${sample}_R2_trim.fq.gz" | gzip > "$BUILD_ROOT/revcomp/${sample}_trim_revcomp.fq.gz"
done
