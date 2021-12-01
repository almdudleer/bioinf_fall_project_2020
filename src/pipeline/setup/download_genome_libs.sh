#!/bin/bash
set -e

[ -n "$BUILD_ROOT" ] || { echo "BUILD_ROOT variable is not set"; exit 1; }

mkdir -p "$BUILD_ROOT/downloads"

if ! [ -f "$BUILD_ROOT/downloads/mature.fa.gz" ] ; then
  wget https://www.mirbase.org/ftp/CURRENT/mature.fa.gz -P "$BUILD_ROOT/downloads"
else
  echo "mature.fa.gz already exists"
fi

if ! [ -f "$BUILD_ROOT/downloads/hsa.gff3" ] ; then
  wget https://www.mirbase.org/ftp/CURRENT/genomes/hsa.gff3 -P "$BUILD_ROOT/downloads"
else
  echo "hsa.gff3 already exists"
fi

if ! [ -f "$BUILD_ROOT/downloads/GRCh38_no_alt.zip" ] ; then
  wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/GRCh38_no_alt.zip -P "$BUILD_ROOT/downloads"
else
  echo "GRCh38_no_alt.zip already exists"
fi

echo "All the necessary genome libraries are downloaded"