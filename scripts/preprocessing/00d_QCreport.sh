#!/usr/bin/env bash

# Containing of each script in the workflow

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PROJECT_ROOT}/config/paths.sh"

# GSE Accession Values: 
# GEO: GSE111657, GSE111658, GSE111659, GSE111661
# PRJNA: PRJNA437760

# fastqc report
parallel -j 3 --eta --bar \
  '
  fastqc -t 3 "{}" && echo "{/.} complete"
  ' \
::: ./*sample*.fastq.gz

