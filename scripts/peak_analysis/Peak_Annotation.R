# Peak Annotation for ChIP using ChIPseeker


# Load libraries
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

# Load in peak files
sample_files <- list.files(path = "./data/macs2/narrowPeak/", full.names = T)
sample_names <- basename(sample_files)
sample_names <- sub("_peaks\\.narrowPeak$", "", sample_names)

for(r in 1:length(sample_files)){
  obj <- ChIPpeakAnno::toGRanges(sample_files[r], format="narrowPeak", header=FALSE)  
  assign(sample_names[r], obj)
}

# Load genome reference for annotation
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Peak Annotation
annot_WT1 <- annotatePeak(WT_H3K27ac_ChIPseq_REP1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

# View the result
annot_WT1@anno %>% 
  data.frame() %>% View()

# Piechart
plotAnnoPie(annot_WT1)

# Barplot
plotAnnoBar(annot_WT1)

# UpSet plot
upsetplot(annot_WT1)

# TSS distance plot
plotDistToTSS(annot_WT1)


## Visualization for multiple samples ##

# Create a list of GRanges objects
samples_list <- list(
  WT1 = WT_H3K27ac_ChIPseq_REP1,
  WT2 = WT_H3K27ac_ChIPseq_REP2,
  WT3 = WT_H3K27ac_ChIPseq_REP3,
  cKO1 = cKO_H3K27ac_ChIPseq_REP1,
  cKO2 = cKO_H3K27ac_ChIPseq_REP2,
  cKO3 = cKO_H3K27ac_ChIPseq_REP3)

# Annotate each sample
peakAnnoList <- lapply(samples_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), annoDb="org.Mm.eg.db", verbose=FALSE)

# Annotate each sample
peakAnnoList <- lapply(samples_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), annoDb="org.Mm.eg.db", verbose=FALSE)

# Plot barplot for all samples
plotAnnoBar(peakAnnoList)

# Plot distance to TSS for all samples
plotDistToTSS(peakAnnoList)


## Enrichment around TSS ##

# Get promoters
promoter <- getPromoters(TxDb = txdb, upstream = 2000, downstream = 2000)
promoter

# Create tag matrix (computationally expensive deepTools profile analagous)
tagMatrix_wt1 <- getTagMatrix(WT_H3K27ac_ChIPseq_REP1, windows = promoter)

# Draw a profile plot for WT rep1
plotAvgProf(tagMatrix_wt1, xlim=c(-2000, 2000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency",
            conf = 0.95, resample = 1000)

# Plot heatmap - this may not plot as it requires a lot of memory!
tagHeatmap(tagMatrix_wt1)

## Plot Profiles for all samples (Requires larger RAM and CPU) ##

# Create a tagMatrix for each sample
#tagMatrixList <- lapply(samples_list, getTagMatrix, windows=promoter)

# Plot profile plots, multiple lines in a single plot
#plotAvgProf(tagMatrixList, xlim=c(-2000, 2000))

# Plot profile plots, faceted by row
#plotAvgProf(tagMatrixList, xlim=c(-2000, 2000), conf=0.95,resample=500, facet="row")

# Plot multiple heatmaps (huge memory requirement)
#tagHeatmap(tagMatrixList)









