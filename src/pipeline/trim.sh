#!/bin/bash
set -e

[ -n "$BUILD_ROOT" ] || { echo "BUILD_ROOT is not set"; exit 1; }
[ -n "$DATA_ROOT" ] || { echo "DATA_ROOT is not set"; exit 1; }

mkdir -p "$BUILD_ROOT/trimmed"

for filename in "$DATA_ROOT"/fastq/*.fastq.gz; do
    IFS='_' read -r -a name_parts <<< "$(basename "$filename" .fastq.gz)"
    sample=${name_parts[0]}
    java -jar "$BUILD_ROOT/lib/Trimmomatic-0.39/trimmomatic-0.39.jar" \
      SE "$filename" "$BUILD_ROOT/trimmed/${sample}_R2_trim.fq.gz" \
      ILLUMINACLIP:"$DATA_ROOT"/adapters_rev_compl.fa:2:30:10
done
