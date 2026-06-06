#!/bin/bash

# Create directories
cd ./results/callpeaks
mkdir -p ./visualization/

# Combine Peaks based off of IDR returns
# bedtools intersect \
# -wo -f 0.3 -r \
# -a ./wt_sample1_peaks_filtered.narrowPeak \
# -b ./wt_sample2_peaks_filtered.narrowPeak \
# > wt_combinedpeaks.bed

# generate bigwig files
bamCoverage -b ~/chipseq_workshop/results/bowtie2/wt_sample2_chip_final.bam \
-o ~/chipseq_workshop/results/visualization/bigWig/wt_sample2_chip.bw \
--binSize 20

# or generate bigwig with chip versus input bamCompare (default --operation is in log2 ratio)
bamCompare -b1 ~/chipseq_workshop/results/bowtie2/wt_sample2_chip_final.bam \
-b2 ~/chipseq_workshop/results/bowtie2/wt_sample2_input_final.bam \
-o ~/chipseq_workshop/results/visualization/bigWig/wt_sample2_chip.bw \
--binSize 20


# Cleanup
mv ../results/macs2/*.log ../logs/