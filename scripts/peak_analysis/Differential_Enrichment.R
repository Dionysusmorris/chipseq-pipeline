# Differential Enrichment with Diffbind

library(DiffBind)
library(edgeR)
library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)
library(ChIPpeakAnno)
library(UpSetR)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(ggrepel)
library(ChIPpeakAnno)
library(ggupset)

# Read in samplesheet metadata
samples <- read.csv("data/DiffBind/metadata.csv")
names(samples)

### DO NOT RUN THIS CODE AS IT REQUIRES YOUR BAM FILES TO BE ACCESSIBLE LOCALLY ###

## Read in data files to create DiffBind object 
dbObj <- dba(sampleSheet=samples, scoreCol=5)

# Count reads to create affinity binding matrix
dbObj <- dba.count(dbObj, bParallel=FALSE) # This is the most computationally intensive part

######################

# Load in the existing DiffBind object
dbObj <- readRDS("data/DiffBind/dbObj.rds")

# Explore 
dbObj

# Create a dataframe with total number of reads in our affinity matrix
info <- dba.show(dbObj)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP, peakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID

libsizes

# PCA plot
dba.plotPCA(dbObj, attributes=DBA_CONDITION, label=DBA_ID, score = DBA_SCORE_NORMALIZED, labelSize = 0.6)

# Plot correlation heatmap
dba.plotHeatmap(dbObj, ColAttributes = DBA_TISSUE,
                score = DBA_SCORE_NORMALIZED)


### Differential Binding Affinity Analysis ###

# Set contrasts
dbObj <- dba.contrast(dbObj, categories = DBA_CONDITION)
dbObj

# Identify differentially bound regions (I usually use Deseq2: DBA_DESEQ2)
dbObj <- dba.analyze(dbObj, method = DBA_ALL_METHODS, bGreylist = FALSE, bBlacklist = FALSE)

# Extract summary and compare across different FDRs
de_summary <- dba.show(dbObj, bContrasts = T, th=0.05)
de_summary
de_summary <- dba.show(dbObj, bContrasts = T, th=0.10)
de_summary
de_summary <- dba.show(dbObj, bContrasts = T, th=0.01)
de_summary


# Overlap between edgeR and DESeq2 results to compare
dba.plotVenn(dbObj, contrast = 1, method = DBA_ALL_METHODS)

# Plot PCA using only DE regions
dba.plotPCA(dbObj, contrast=1, method=DBA_DESEQ2, attributes=DBA_CONDITION, label=DBA_ID)

# MA plot
dba.plotMA(dbObj, method=DBA_DESEQ2)

# Volcano plot
dba.plotVolcano(dbObj, contrast = 1)

# Plot heatmap
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
readscores <- dba.plotHeatmap(dbObj, correlations = FALSE,
                              scale="row", colScheme = hmap)

# Extract results
res_deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)

# Write GRanges to rds
saveRDS(res_deseq, file = "all_res_deseq2.rds")

# Write results to a tsv file
out <- as.data.frame(res_deseq)
write.table(out, file="results/cKO_vs_WT_deseq2.txt", sep="\t", quote=F, row.names=F)

# Create bed files for each keeping only significant peaks (p < 0.05)
cKO_enrich <- out %>% 
  filter(FDR < 0.05 & Fold > 0) %>% 
  dplyr::select(seqnames, start, end)
WT_enrich <- out %>% 
  filter(FDR < 0.05 & Fold < 0) %>% 
  dplyr::select(seqnames, start, end)

# Write to file
write.table(cKO_enrich, file="results/cKO_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.table(WT_enrich, file="results/WT_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)


# NOTE: BED files cannot contain headers and so we have added the col.names=F argument
# to address that. Additionally, we took only the first three columns
# from the results (genomic coordinates) to adhere to a minimal BED file format.
