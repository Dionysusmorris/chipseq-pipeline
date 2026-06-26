#!/usr/bin/env bash

# Path variables used to source into the workflow scripts

# Top level script variables
PROJECT_ROOT="."

ORGANISM="mm10"

# Metadata 
META_DIR="${PROJECT_ROOT}/metadata/"

# Fetch SRA fastq
RAW_DATA_DIR="${PROJECT_ROOT}/data/raw"
PROCESSED_DATA_DIR="${PROJECT_ROOT}/data/processed"
SCRATCH="/dev/shm/tmp"

# Reference paths
REF_DIR="${PROJECT_ROOT}/data/reference/${ORGANISM}}"
GENOME="${REF_DIR}/genome/${ORGANISM}.fa"
INDEX="${REF_DIR}/bowtie2/${ORGANISM}"

# Alignment
RESULTS_DIR="${PROJECT_ROOT}/results"
BAM_DIR="${RESULTS_DIR}/bam"
PEAK_DIR="${RESULTS_DIR}/peak_calling"

# Alignment configs
THREADS=8
MAPQ=30

# QC report
QC_DIR="${RESULTS_DIR}/qc"

# peak calling and idr



# peak processing



# generating peak analysis qc metrics





