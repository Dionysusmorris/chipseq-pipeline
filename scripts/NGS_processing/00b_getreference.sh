#!/bin/bash
#

# GSE Accession Values: 
# GEO: GSE111657, GSE111658, GSE111659, GSE111661
# PRJNA: PRJNA437760

# I typically retrieve references via FTP from Ensembl (Alternatives are igenome, GENCODE, or EMBL)
# Ensembl Table of References: https://useast.ensembl.org/info/data/ftp/index.html?
# Ensemble FTPs of nonvertebrates: http://ftp.ensemblgenomes.org/pub/
# Ensemble FTP: https://ftp.ensembl.org/

mkdir -p ./reference_data/
cd ./reference_data/

# fetched reference from igenomes
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz

tar -I "pigz -p 20" -xvf Mus_musculus_UCSC_mm10.tar.gz

# original tar extraction was 21 minutes ..
# using 20 cores and 500mbps I/O speed took 2 minutes

#bowtie2-build <path_to_reference_genome.fa> <prefix_to_name_indexes>

# Filter blacklist regions
# https://github.com/Boyle-Lab/Blacklist/blob/61a04d2c5e49341d76735d485c61f0d1177d08a8/lists/mm10-blacklist.v2.bed.gz