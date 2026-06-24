# R version: 4.6.0
# Bioconductor version: 3.23
# Generated on 2026-06-05


install.packages("BiocManager")
install.packages("remotes")
install.packages("optparse")

remotes::install_github("bcbio/bcbioR")
# CRAN packages

install.packages("BiocManager")
install.packages("tidyverse")
install.packages("pheatmap")
install.packages("UpSetR")
install.packages("RColorBrewer")
install.packages("ggrepel")
install.packages("ggupset")

BiocManager::install("ChIPseeker")
BiocManager::install("ChIPpeakAnno")
BiocManager::install("DiffBind")
BiocManager::install("clusterProfiler")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("IRanges")
BiocManager::install("GenomicRanges")
BiocManager::install("DESeq2")
BiocManager::install("org.Mm.eg.db")