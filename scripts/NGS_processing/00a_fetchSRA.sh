#!/bin/bash
#

# GSE Accession Values: 
# GEO: GSE111657, GSE111658, GSE111659, GSE111661
# PRJNA: PRJNA437760

cd ./raw_data/
rm -rf /dev/shm/tmp
mkdir -p /dev/shm/tmp

parallel -j 3 --delay 2 '
    echo "===== {} =====" &&
    prefetch -p {} -O ./ &&
    fasterq-dump --threads 8 {} --temp /dev/shm/tmp
' :::: "../SRR_samples.txt" && echo "SRA samples fetched"

# Fasterq pull and compression. Here you can merge your techincal replicates that each have unique SRR IDs
# (e.g. Sample A Replicate 1, Replicate 2, Replicate 3)
# Rework this snippet. It is incorrect.
#
# parallel -j 3 --delay 2 '
#     echo "===== {} =====" &&
#     prefetch -p {} -O ./ &&
#     fasterq-dump --threads 8 {} --temp /dev/shm/tmp &&
#     wait
#     pigz -p 15 -c {} > $some_name
# ' :::: "../SRR_samples.txt" && echo "SRA samples fetched"


# echo merging and gzipping fastqs
pigz -p 15 -c SRR6823762.fastq SRR6823763.fastq > wt_sample1_chip.fastq.gz
pigz -p 15 -c SRR6823764.fastq SRR6823765.fastq > wt_sample1_input.fastq.gz
pigz -p 15 -c SRR6823766.fastq SRR6823767.fastq > wt_sample2_chip.fastq.gz
pigz -p 15 -c SRR6823768.fastq SRR6823769.fastq > wt_sample2_input.fastq.gz

pigz -p 15 -c SRR6823770.fastq > ko_sample1_chip.fastq.gz
pigz -p 15 -c SRR6823771.fastq > ko_sample1_input.fastq.gz
pigz -p 15 -c SRR6823772.fastq > ko_sample2_chip.fastq.gz
pigz -p 15 -c SRR6823773.fastq > ko_sample2_input.fastq.gz

rm ./*.fastq
rm -rf /dev/shm/tmp
mkdir -p /dev/shm/tmp

cd ../
echo complete

# Cleanups
# rm -rf ./raw_data/SRR*
# rm -rf ./raw_data/*.fastq

