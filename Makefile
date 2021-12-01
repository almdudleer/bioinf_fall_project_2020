.PHONY: all clean clean_all setup_tools qc
SHELL:=/bin/bash

export SRC_ROOT := src
export BUILD_ROOT := build
export DATA_ROOT := data
export SAMPLES := N3 N18 N20 N22 N27 N31 N37 N47 N48 N51 N58 N70 N74

GENOME_TAGGED=$(patsubst %, build/genome/alignment/%_aligned_genome_tagged.bam, $(SAMPLES))
GENOME_RENAMED=$(patsubst %, build/genome/alignment/%_aligned_genome_renamed.bam, $(SAMPLES))
GENOME_SORTED=$(patsubst %, build/genome/alignment/%_aligned_genome.bam, $(SAMPLES))
GENOME_ALIGNED=$(patsubst %, build/genome/alignment/%_aligned_genome.sam, $(SAMPLES))

MATURE_UNALIGNED=$(patsubst %, build/mature/alignment/%_unaligned_maturemiRNA.fq, $(SAMPLES))
MATURE_SORTED=$(patsubst %, build/mature/alignment/%_aligned_maturemiRNA.bam, $(SAMPLES))
MATURE_ALIGNED=$(patsubst %, build/mature/alignment/%_aligned_maturemiRNA.sam, $(SAMPLES))


# Phony targets
all: build/out/genome_counts.tsv build/out/mirna_counts.tsv

qc: build/out/genome_counts.tsv build/out/mirna_counts.tsv
	cd build && multiqc .

# Genome
build/out/genome_counts.tsv: build/genome/counts
	./src/pipeline/genome/combine_counts.sh

build/genome/counts: $(GENOME_TAGGED)
	./src/pipeline/genome/get_counts.sh

$(GENOME_TAGGED): $(GENOME_RENAMED)
	./src/pipeline/genome/tag.sh

$(GENOME_RENAMED): $(GENOME_SORTED)
	./src/pipeline/genome/reheader.sh

$(GENOME_SORTED): $(GENOME_ALIGNED)
	./src/pipeline/genome/sort.sh

$(GENOME_ALIGNED): $(MATURE_UNALIGNED) build/genome/index
	./src/pipeline/genome/align.sh

build/genome/index: build/downloads/GRCh38_no_alt.zip
	./src/pipeline/genome/index.sh

# Mature miRNA alignment
build/out/mirna_counts.tsv: build/mature/counts
	./src/pipeline/mature/combine_counts.sh

build/mature/counts: $(MATURE_SORTED)
	./src/pipeline/mature/get_counts.sh

$(MATURE_SORTED): $(MATURE_ALIGNED)
	./src/pipeline/mature/sort.sh

$(MATURE_ALIGNED) $(MATURE_UNALIGNED): build/short_long build/mature/index
	./src/pipeline/mature/align.sh

build/mature/index: build/downloads/mature.fa.gz
	./src/pipeline/mature/index.sh

# Common part
build/short_long: build/revcomp
	./src/pipeline/short_long_split.sh

build/revcomp: build/trimmed
	./src/pipeline/revcomp.sh

build/trimmed: data/fastq build/lib/Trimmomatic-0.39/trimmomatic-0.39.jar
	./src/pipeline/trim.sh

build/lib/Trimmomatic-0.39/trimmomatic-0.39.jar:
	./src/pipeline/setup/setup_tools.sh

build/downloads/GRCh38_no_alt.zip build/downloads/mature.fa.gz:
	./src/pipeline/setup/download_genome_libs.sh

# Utils
setup_tools:
	ret=0
	src/pipeline/setup/check_tools.sh || ret=$?
	if [ $ret -eq 2 ] ; then
	  echo "Trying to install missing utils automatically"
	  src/pipeline/setup/setup_tools.sh
	fi
	src/pipeline/setup/check_data.sh

clean:
	find build ! -name 'build' ! -path 'build/downloads*' -exec rm -rf {} +

clean_all:
	rm -rf build
