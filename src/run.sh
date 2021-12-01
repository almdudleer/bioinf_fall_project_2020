#!/bin/bash
set -e

SRC_ROOT=.
BUILD_ROOT=../build
DATA_ROOT=../data
SAMPLES="N3 N18 N20 N22 N27 N31 N37 N47 N48 N51 N58 N70 N74"
export SAMPLES DATA_ROOT SRC_ROOT BUILD_ROOT


# Setup
## Perform checks
ret=0
"$SRC_ROOT/pipeline/setup/check_tools.sh" || ret=$?
if [ $ret -eq 2 ] ; then
  echo "Trying to install missing utils automatically"
  "$SRC_ROOT/pipeline/setup/setup_tools.sh"
fi
"$SRC_ROOT/pipeline/setup/check_data.sh"

## Download mature miRNA and human genome references
"$SRC_ROOT/pipeline/setup/download_genome_libs.sh"


# Run the pipeline
## Common part
"$SRC_ROOT/pipeline/trim.sh"
"$SRC_ROOT/pipeline/revcomp.sh"
"$SRC_ROOT/pipeline/short_long_split.sh"

## Mature miRNA alignment part
"$SRC_ROOT/pipeline/mature/index.sh"
"$SRC_ROOT/pipeline/mature/align.sh"
"$SRC_ROOT/pipeline/mature/sort.sh"
"$SRC_ROOT/pipeline/mature/get_counts.sh"
"$SRC_ROOT/pipeline/mature/combine_counts.sh"


## Genome alignment part
"$SRC_ROOT/pipeline/genome/index.sh"
"$SRC_ROOT/pipeline/genome/align.sh"
"$SRC_ROOT/pipeline/genome/sort.sh"
"$SRC_ROOT/pipeline/genome/reheader.sh"
"$SRC_ROOT/pipeline/genome/tag.sh"
"$SRC_ROOT/pipeline/genome/get_counts.sh"
"$SRC_ROOT/pipeline/genome/combine_counts.sh"
