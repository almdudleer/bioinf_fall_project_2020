#!/bin/bash
set -e

SRC_ROOT=.
BUILD_ROOT=../build
DATA_ROOT=../data
SAMPLES="N3 N18 N20 N22 N27 N31 N37 N47 N48 N51 N58 N70 N74"
export SAMPLES DATA_ROOT SRC_ROOT BUILD_ROOT


# Setup
## Perform checks
tools_installed=0
for util in {wget,zcat,unzip,gzip,java,R,python3,conda}; do
  which "$util" 1>/dev/null || { echo "You need to install $util manually"; tools_installed=1; }
done
[ $tools_installed -eq 0 ] || { echo "Install the utils listed above to continue"; exit 1; }

for util in {seqkit,samtools,bedtools,bowtie,bowtie-build,fastqc,multiqc}; do
  which "$util" 1>/dev/null || { echo "$util is missing" ; tools_installed=1; }
done

if ! [ -e "$BUILD_ROOT/lib/Trimmomatic-0.39/trimmomatic-0.39.jar" ] ; then
  echo "Trimmomatic is missing"
  tools_installed=1
fi

if [ $tools_installed -eq 1 ] ; then
  echo "Trying to install missing utils automatically"
  conda install -c bioconda --yes --file "$SRC_ROOT/requirements.txt"
  echo "Installing trimmomatic"
  mkdir -p "$BUILD_ROOT/downloads"
  wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip -P "$BUILD_ROOT/downloads"
  mkdir -p "$BUILD_ROOT/lib"
  unzip "$BUILD_ROOT/downloads/Trimmomatic-0.39.zip" -d "$BUILD_ROOT/lib"
else
  echo "All the required utils are installed"
fi

all_files_present=0
required_files=("adapters_rev_compl.fa" "ncbitochrom.txt" "mirnas.txt" "hsa-genome-miRBase22v-onlymiRNAs-convforTagBAM.bed")
IFS=' ' read -r -a samples_array <<< "$SAMPLES"

for i in "${!samples_array[@]}" ; do
  index=$(( $i + 1 ))
  required_files+=( fastq/"${samples_array[$i]}"_S"$index"_L001_R2_001.fastq.gz )
done

for file in "${required_files[@]}" ; do
  if ! [ -e "$DATA_ROOT/$file" ] ; then
    echo "$file is missing"
    all_files_present=1
  fi
done

if [ $all_files_present -ne 0 ] ; then
  echo "Put all the missing files into the data directory"
  exit 1
else
  echo "All the required files are present"
fi

## Download mature miRNA and human genome references
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

# Common part
## Trim
mkdir -p "$BUILD_ROOT/trimmed"

for filename in "$DATA_ROOT"/fastq/*.fastq.gz; do
    IFS='_' read -r -a name_parts <<< "$(basename "$filename" .fastq.gz)"
    sample=${name_parts[0]}
    java -jar "$BUILD_ROOT/lib/Trimmomatic-0.39/trimmomatic-0.39.jar" \
      SE "$filename" "$BUILD_ROOT/trimmed/${sample}_R2_trim.fq.gz" \
      ILLUMINACLIP:"$DATA_ROOT"/adapters_rev_compl.fa:2:30:10
done

## Reverse complement
mkdir -p "$BUILD_ROOT/revcomp"

for sample in $SAMPLES ; do
  seqkit seq -r -p "$BUILD_ROOT/trimmed/${sample}_R2_trim.fq.gz" | gzip > "$BUILD_ROOT/revcomp/${sample}_trim_revcomp.fq.gz"
done

## Short long split
mkdir -p "$BUILD_ROOT/short_long"

for sample in $SAMPLES ; do
  filename="$BUILD_ROOT/revcomp/${sample}_trim_revcomp.fq.gz"
  seqkit seq "$filename" -M 31 > "$BUILD_ROOT/short_long/${sample}_short.fq"
  seqkit seq "$filename" -m 32 > "$BUILD_ROOT/short_long/${sample}_long.fq"
  seqkit seq "$BUILD_ROOT/short_long/${sample}_short.fq" -m 15 > "$BUILD_ROOT/short_long/${sample}_short_15-31.fq"
done

# Mature miRNA reference part
## Reference indexing
mkdir -p "$BUILD_ROOT/mature/index"

zcat < "$BUILD_ROOT/downloads/mature.fa.gz" | seqkit grep -r -p ^hsa | seqkit seq --rna2dna > "$BUILD_ROOT/mature/mature-hsa.fa"
bowtie-build "$BUILD_ROOT/mature/mature-hsa.fa" "$BUILD_ROOT/mature/index/mature"

## Alignment
mkdir -p "$BUILD_ROOT/mature/alignment"

for sample in $SAMPLES ; do
    filename="$BUILD_ROOT/short_long/${sample}_short_15-31.fq"
    bowtie -n 0 -l 32 --norc --best --strata -m 1 "$BUILD_ROOT/mature/index/mature" "$filename" \
     --un "$BUILD_ROOT/mature/alignment/${sample}_unaligned_maturemiRNA.fq" \
     -S "$BUILD_ROOT/mature/alignment/${sample}_aligned_maturemiRNA.sam" \
     2> "$BUILD_ROOT/mature/alignment/${sample}_maturemiRNA.log"
done

## Sorting and indexing
for sample in $SAMPLES ; do
    samtools sort "$BUILD_ROOT/mature/alignment/${sample}_aligned_maturemiRNA.sam" > "$BUILD_ROOT/mature/alignment/${sample}_aligned_maturemiRNA.bam"
    samtools index "$BUILD_ROOT/mature/alignment/${sample}_aligned_maturemiRNA.bam"
done

## Getting read counts
mkdir -p "$BUILD_ROOT/mature/counts"

for sample in $SAMPLES ; do
    samtools idxstats "$BUILD_ROOT/mature/alignment/${sample}_aligned_maturemiRNA.bam" | cut -f1,3 - | \
     sed -e "1s/^/miRNA\t${sample}-miRNAcount\n/" > "$BUILD_ROOT/mature/counts/${sample}_mature_counts.tsv"
done

## Combine all sample counts with an R script
mkdir -p "$BUILD_ROOT/out"

"$SRC_ROOT/pipeline/utils/combine_counts.R" "$BUILD_ROOT/mature/counts" "*.tsv" "$BUILD_ROOT/out/mirna_counts.tsv"


# Genome reference part

## Unpacking the pre-built index
mkdir -p "$BUILD_ROOT/genome/index"
unzip "$BUILD_ROOT/downloads/GRCh38_no_alt.zip" -d "$BUILD_ROOT/genome/index"

## Alignment
mkdir -p "$BUILD_ROOT/genome/alignment"

for sample in $SAMPLES ; do
  filename="$BUILD_ROOT/mature/alignment/${sample}_unaligned_maturemiRNA.fq"
  bowtie -n 1 -l 32 --norc --best --strata -m 1 \
   "$BUILD_ROOT/genome/index/GCA_000001405.15_GRCh38_no_alt_analysis_set" "$filename" \
   --al "$BUILD_ROOT/genome/alignment/${sample}_aligned_genome.fastq" \
   -S "$BUILD_ROOT/genome/alignment/${sample}_aligned_genome.sam" \
   2> "$BUILD_ROOT/genome/alignment/${sample}_aligned_genomeaftermiRNA.log"
done

## Sorting
for sample in $SAMPLES ; do
  filename="$BUILD_ROOT/genome/alignment/${sample}_aligned_genome.sam"
  samtools sort "$filename" > "$BUILD_ROOT/genome/alignment/${sample}_aligned_genome.bam"
  samtools index "$BUILD_ROOT/genome/alignment/${sample}_aligned_genome.bam"
done

## Changing header chromosome format
for sample in $SAMPLES ; do
    header=$(samtools view -H "$BUILD_ROOT/genome/alignment/${sample}_aligned_genome.bam")
    map="$DATA_ROOT/ncbitochrom.txt"
    while IFS=$'\t' read -r ucsc_name ncbi_name; do
        header=${header//$ncbi_name/$ucsc_name/}
    done < "$map"
    echo "$header" > "$BUILD_ROOT/genome/alignment/${sample}_tmp_header.txt"
    samtools reheader "$BUILD_ROOT/genome/alignment/${sample}_tmp_header.txt" \
      "$BUILD_ROOT/genome/alignment/${sample}_aligned_genome.bam" \
      > "$BUILD_ROOT/genome/alignment/${sample}_aligned_genome_renamed.bam"
    rm "$BUILD_ROOT/genome/alignment/${sample}_tmp_header.txt"
done

## Tagging
for sample in $SAMPLES ; do
  bedtools tag -i "$BUILD_ROOT/genome/alignment/${sample}_aligned_genome_renamed.bam" \
    -files "$DATA_ROOT/hsa-genome-miRBase22v-onlymiRNAs-convforTagBAM.bed" -names -tag XQ \
    > "$BUILD_ROOT/genome/alignment/${sample}_aligned_genome_tagged.bam"
done

## Getting read counts
mkdir -p "$BUILD_ROOT/genome/counts"

for sample in $SAMPLES ; do
  filename="$BUILD_ROOT/genome/alignment/${sample}_aligned_genome_tagged.bam"
  printf "\nProcessing sample: %s\n" "$sample"
  samtools view "$filename" | \
    "$SRC_ROOT/pipeline/utils/get_counts.py" -m "$DATA_ROOT/mirnas.txt" -o "$BUILD_ROOT/genome/counts/${sample}_genome_based_counts.txt"
  sed "N;s/\n/\t/g" "$BUILD_ROOT/genome/counts/${sample}_genome_based_counts.txt" | \
   sort -k1 | uniq | sed "1s/^/miRNA\t${sample}\n/" > "$BUILD_ROOT/genome/counts/${sample}_taggedBAMcounts.txt"
done

## Combining the counts
mkdir -p "$BUILD_ROOT/out"

"$SRC_ROOT/pipeline/utils/combine_counts.R" "$BUILD_ROOT/genome/counts" "*_taggedBAMcounts.txt" "$BUILD_ROOT/out/genome_counts.tsv"
