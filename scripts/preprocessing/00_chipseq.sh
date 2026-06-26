#!/usr/bin/env bash

# I prefer executing this script at the top level of your project directory
# bash WORKING_PROJECT_DIRECTORY/scripts/preprocessing/00_chipseq.sh

# Establish paths for project if not already existing
mkdir -p \
    ${PROJECT_ROOT}/data/{raw,processed/trimmed,reference/{genome,annotation,blacklist,bowtie2,bwa}} \
    results/{qc/{fastqc,multiqc},bam,peak_calling/{idr},counts,differential,figures,tables} \
    {config,metadata,scripts,envs,docs,logs,nextflow,tmp}

# Absolute path to this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Project root (two directories above this script)
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
# Source project paths as variables for simple recall
source "${PROJECT_ROOT}/config/paths.sh"

# Export of variables for end-to-end workflow
export GENOME
export INDEX
export ORGANISM




#### Run each step in the workflow #####

# retrieve and parse metadata (To be completed)

# fetch SRA raw fastqs
bash "${PREPROCESS_DIR}/00a_fetchSRA.sh"

# fastq trimming script (To be completed)

# get reference genome
bash "${PREPROCESS_DIR}/00b_getreference.sh"

# run alignment
bash "${PREPROCESS_DIR}/00c_runbowtie2.sh"

# generate alignment QC Report 
bash "${PREPROCESS_DIR}/00d_QCreport.sh"

# call peaks
bash "${PREPROCESS_DIR}/00e_callpeaks.sh"

# generate idr report
bash "${PREPROCESS_DIR}/00f_idr.sh"

# merge concording peaks
bash "${PREPROCESS_DIR}/00g_peakmerge.sh"

# generate peak calling qc report
bash "${PREPROCESS_DIR}/00h_generate_qcreport.sh"

# generate post peakcalling metrics for peak analysis qc
bash "${PREPROCESS_DIR}/00z_calculate_frip.sh"

