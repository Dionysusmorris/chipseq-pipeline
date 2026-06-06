#!/bin/bash

# Create directories
cd ./results/
mkdir -p ./callpeaks/

# Call peaks with MACS3

parallel -j 4 --eta --bar '
    set -euo pipefail

    chip="{}"
    base=$(basename "$chip")
    name=${base%_chip.sorted.bam}

    input="./bam/sorted/${name}_input.sorted.bam"

    macs3 callpeak \
    -t "$chip" \
    -c "$input" \
    --broad \
    --broad-cutoff 0.1 \
    -f BAM \
    -g mm \
    -n "$name" \
    --outdir ./callpeaks/ \
    2> ./callpeaks/${name}.log

    echo "${name} peakset generated"
' ::: ./bam/sorted/*_chip.sorted.bam


# Cleanup

mv ../results/callpeaks/*.log ../logs/

