#!/usr/bin/env bash

# Containing of each script in the workflow

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PROJECT_ROOT}/config/paths.sh"

# GSE Accession Values: 
# GEO: GSE111657, GSE111658, GSE111659, GSE111661
# PRJNA: PRJNA437760

cd "${RAW_DATA_DIR}"
rm -rf "${SCRATCH}"
mkdir -p "${SCRATCH}"

parallel \
    --env RAW_DATA_DIR,SCRATCH,THREADS \
    -j 3 \
    --delay 2 \
    '
    echo "===== {} =====" &&
    prefetch -p {} -O "${RAW_DATA_DIR}" &&
    fasterq-dump \
    --threads "${THREADS}" {} \
    --temp "${SCRATCH}" \
    --outdir "${RAW_DATA_DIR}"
    ' :::: "${META_DIR}/SRR_samples.txt" && echo "SRA samples fetched"

# Fasterq pull and compression. Here you can merge your techincal replicates that each have unique SRR IDs
# (e.g. Sample A Replicate 1, Replicate 2, Replicate 3)
# Rework this snippet. It is incorrect.
#
# parallel -j 3 --delay 2 '
#     echo "===== {} =====" &&
#     prefetch -p {} -O "${RAW_DATA_DIR}" &&
#     fasterq-dump --threads 8 {} --temp "${SCRATCH}"" &&
#     wait
#     pigz -p 15 -c {} > $some_name
# ' :::: "${META_DIR}/SRR_samples.txt" && echo "SRA samples fetched"


# echo merging and gzipping fastqs
pigz -p 15 -c SRR6823762.fastq SRR6823763.fastq > wt_sample1_chip.fastq.gz
pigz -p 15 -c SRR6823764.fastq SRR6823765.fastq > wt_sample1_input.fastq.gz
pigz -p 15 -c SRR6823766.fastq SRR6823767.fastq > wt_sample2_chip.fastq.gz
pigz -p 15 -c SRR6823768.fastq SRR6823769.fastq > wt_sample2_input.fastq.gz

pigz -p 15 -c SRR6823770.fastq > ko_sample1_chip.fastq.gz
pigz -p 15 -c SRR6823771.fastq > ko_sample1_input.fastq.gz
pigz -p 15 -c SRR6823772.fastq > ko_sample2_chip.fastq.gz
pigz -p 15 -c SRR6823773.fastq > ko_sample2_input.fastq.gz

rm "${RAW_DATA_DIR}/*.fastq"
rm -rf "${SCRATCH}"
mkdir -p "${SCRATCH}"

cd "${PROJECT_DIR}"
echo complete

# Cleanups
# rm -rf "${RAW_DATA_DIR}/SRR*"
# rm -rf "${RAW_DATA_DIR}/*.fastq"

