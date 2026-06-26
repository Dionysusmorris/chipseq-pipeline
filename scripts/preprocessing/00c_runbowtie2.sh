#!/usr/bin/env bash

# Containing of each script in the workflow

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PROJECT_ROOT}/config/paths.sh"

# GSE Accession Values: 
# GEO: GSE111657, GSE111658, GSE111659, GSE111661
# PRJNA: PRJNA437760


# ALIGN, CONVERT TO BINARY, FILTER, SORT

parallel -j 3 --joblog ../logs/alignment.log --eta --bar '
  set -euo pipefail
  f="{}";
  base=$(basename "$f");
  name=${base%.fastq.gz};

  time bowtie2 --threads 8 -q --local \
    -x ../reference_data/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome \
    -U "$f" \
    2> ../logs/${name}.log \
  | samtools view -@ 2 -b -q 30 -F 260 - \
  | bedtools intersect -v -abam - -b ../reference_data/mm10-blacklist.v2.bed \
  | samtools sort -@ 4 -m 4G -T ../tmp/${name}.tmp -o ./bam/sorted/${name}.sorted.bam - \
  && samtools index -@ 2 ./bam/sorted/${name}.sorted.bam;

  echo "${name} complete"
' ::: ../raw_data/*.fastq.gz

# Duplicate Mapped Reads are not filtered here

# Cleanup
rm -r ../tmp
cd ../
