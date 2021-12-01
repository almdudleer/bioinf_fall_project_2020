#!/bin/bash
set -e

[ -n "$SAMPLES" ] || { echo "SAMPLES is not set"; exit 1; }
[ -n "$DATA_ROOT" ] || { echo "DATA_ROOT is not set"; exit 1; }

exit_code=0

required_files=("adapters_rev_compl.fa" "ncbitochrom.txt" "mirnas.txt" "hsa-genome-miRBase22v-onlymiRNAs-convforTagBAM.bed")
IFS=' ' read -r -a samples_array <<< "$SAMPLES"

for i in "${!samples_array[@]}" ; do
  index=$(( $i + 1 ))
  required_files+=( fastq/"${samples_array[$i]}"_S"$index"_L001_R2_001.fastq.gz )
done

for file in "${required_files[@]}" ; do
  if ! [ -e "$DATA_ROOT/$file" ] ; then
    echo "$file is missing"
    exit_code=1
  fi
done

if [ $exit_code -ne 0 ] ; then
  echo "Put all the missing files into the data directory"
else
  echo "All the required files are present"
fi

exit $exit_code
