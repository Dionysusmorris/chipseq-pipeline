#!/bin/bash

# Create directories
cd ./results
mkdir -p ./figures/ ./bed/

# Combine Peaks based off of IDR returns
# bedtools intersect \
# -wo -f 0.3 -r \
# -a ./callpeaks/wt_sample1_peaks_filtered.narrowPeak \
# -b ./callpeaks/wt_sample2_peaks_filtered.narrowPeak \
# > ./bed/wt_combinedpeaks.bed

# generate bigwig files
bamCoverage -b cd ./bam/sorted/wt_sample2_chip_final.bam \
-o ./results/figures/bigWig/wt_sample2_chip.bw \
--binSize 20

# or generate bigwig with chip versus input bamCompare (default --operation is in log2 ratio)
bamCompare -b1 ./bam/sorted/wt_sample2_chip_final.bam \
-b2 ./bam/sorted/wt_sample2_input_final.bam \
-o ./results/figures/bigWig/wt_sample2_chip.bw \
--binSize 20


# Cleanup
mv ../results/macs2/*.log ../logs/