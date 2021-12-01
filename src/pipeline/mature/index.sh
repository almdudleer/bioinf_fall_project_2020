#!/bin/bash
set -e

[ -n "$BUILD_ROOT" ] || { echo "BUILD_ROOT is not set"; exit 1; }

mkdir -p "$BUILD_ROOT/mature/index"

zcat < "$BUILD_ROOT/downloads/mature.fa.gz" | seqkit grep -r -p ^hsa | seqkit seq --rna2dna > "$BUILD_ROOT/mature/mature-hsa.fa"
bowtie-build "$BUILD_ROOT/mature/mature-hsa.fa" "$BUILD_ROOT/mature/index/mature"
