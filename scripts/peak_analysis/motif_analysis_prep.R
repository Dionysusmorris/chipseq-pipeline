## Preparing data for motif analysis

# Load libraries
library(tidyverse)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# Get merged peaks for WT
wt_consensus <- olaps_wt$mergedPeaks

# Filter to keep only regions 1000bp or smaller
wt_consensus <- wt_consensus[width(wt_consensus) <= 1000]

# Remove non standard sequences as it will cause errors when retrieving the FASTA sequences from the reference (i.e. chr1_GL456211_random)
wt_consensus <- keepStandardChromosomes(wt_consensus, pruning.mode="coarse")

# Annotate consensus reiogns 
annot_wt_consensus <- annotatePeak(wt_consensus, tssRegion = c(-3000, 3000),
                                   TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                                   annoDb = "org.Mm.eg.db")

# GRanges to dataframe
df <- annot_wt_consensus@anno %>% data.frame()

# Find regions annotated as "Promoter" and remove them
promoter_regions <- grep("Promoter", df$annotation)
non_promoter_df <- df[-(promoter_regions), ]

# Write minimal BED file
write_tsv(non_promoter_df[, 1:3], file = "results/non_promoter_wt.bed", 
          col_names = F, quote = "none")



