#!/bin/bash
exit_code=0
for util in {wget,zcat,unzip,gzip,java,R,python3,conda}; do
  which "$util" 1>/dev/null || { echo "You need to install $util manually"; exit_code=1; }
done
[ $exit_code -eq 0 ] || { echo "Install the utils listed above to continue"; exit $exit_code; }

for util in {seqkit,samtools,bedtools,bowtie,bowtie-build,fastqc,multiqc}; do
  which "$util" 1>/dev/null || { echo "$util is missing" ; exit_code=2; }
done

if ! [ -e "$BUILD_ROOT/lib/Trimmomatic-0.39/trimmomatic-0.39.jar" ] ; then
  echo "Trimmomatic is missing"
  exit_code=2
fi

if [ $exit_code -eq 2 ] ; then
  echo "Missing utils can be installed automatically with setup_tools.sh"
else
  echo "All the required utils are installed"
fi

exit $exit_code
