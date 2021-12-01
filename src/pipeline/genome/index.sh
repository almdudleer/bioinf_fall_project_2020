#!/bin/bash
set -e

[ -n "$BUILD_ROOT" ] || { echo "BUILD_ROOT is not set"; exit 1; }

mkdir -p "$BUILD_ROOT/genome/index"

# There will be no actual indexing, since the downloaded genome is already indexed
unzip "$BUILD_ROOT/downloads/GRCh38_no_alt.zip" -d "$BUILD_ROOT/genome/index"