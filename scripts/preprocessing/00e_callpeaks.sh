#!/usr/bin/env bash

# Containing of each script in the workflow

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PROJECT_ROOT}/config/paths.sh"

# GSE Accession Values: 
# GEO: GSE111657, GSE111658, GSE111659, GSE111661
# PRJNA: PRJNA437760
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
