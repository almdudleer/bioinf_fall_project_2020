#!/bin/bash
set -e

[ -n "$SRC_ROOT" ] || { echo "SRC_ROOT is not set"; exit 1; }
[ -n "$BUILD_ROOT" ] || { echo "BUILD_ROOT variable is not set"; exit 1; }

conda install -c bioconda --yes --file "$SRC_ROOT/requirements.txt"

echo "Installing trimmomatic"

mkdir -p "$BUILD_ROOT/downloads"
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip -P "$BUILD_ROOT/downloads"

mkdir -p "$BUILD_ROOT/lib"
unzip "$BUILD_ROOT/downloads/Trimmomatic-0.39.zip" -d "$BUILD_ROOT/lib"
