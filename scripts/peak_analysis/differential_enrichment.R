# Differential Enrichment with Diffbind

# example:
# Rscript differential_enrichment.R \
# --meta_csv meta/diffbind_metadata.csv \
# --output_dir results/differential_enrichment

# DiffBind requires you to create a metadata file:

# Load Libraries
library(DiffBind)
library(tidyverse)
library(optparse)


#### INPUT ARGUMENT OPTIONS FOR COMMANDLINE #####

option_list <- list(
  make_option("--meta_csv", type = "character"),
  make_option("--output_dir", type = "character")
  )

opt <- parse_args(OptionParser(option_list = option_list))
output_dir <- opt$output_dir
meta_csv <- opt$meta_csv
samples <- names(meta_csv)


### REQUIRES YOUR BAM FILES TO BE ACCESSIBLE ###

## Read in BAM files to create DiffBind object
dbObj <- dba(sampleSheet=meta_csv, scoreCol=5)

# Explore
dbObj

# Count reads to create affinity binding matrix

#register(SnowParam(4))
system.time({
  dbObj <- dba.count(dbObj, bParallel = FALSE) # took 15 minutes
}) # This is the most computationally intensive part and may want to scale to CPU and RAM power

# Explore
dbObj

# Create a dataframe with total number of reads in our affinity matrix
info <- dba.show(dbObj)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP, peakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID

libsizes


#############################


# PCA plot

pdf(
  file.path(output_dir, "DiffBind_PCA.pdf"),
  width = 8,
  height = 6
)

dba.plotPCA(dbObj, 
            attributes=DBA_CONDITION, 
            label=DBA_ID, 
            score = DBA_SCORE_NORMALIZED, 
            labelSize = 0.6
            )

dev.off()


# Plot correlation heatmap

pdf(
  file.path(output_dir, "DiffBind_sample_correlation_heatmap.pdf"),
  width = 8,
  height = 8
)

dba.plotHeatmap(
  dbObj,
  ColAttributes = DBA_TISSUE,
  score = DBA_SCORE_NORMALIZED
)

dev.off()

message(
  "Saved differential peaks PCA & heatmap: ",
  file.path(output_dir, "Sample_correlation_heatmap.pdf"),
  file.path(output_dir, "DiffBind_sample_correlation_heatmap.pdf")
)


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



##### PLOTS #####

# Overlap between edgeR and DESeq2 results to compare
pdf(
  file.path(output_dir, "Bindsite_overlaps_venndiagram.pdf"),
  width = 8,
  height = 8
)

dba.plotVenn(dbObj, 
             contrast = 1,
             method = DBA_ALL_METHODS
             )

dev.off()

# Plot PCA using only DE regions

pdf(
  file.path(output_dir, "Differential_enriched_PCA.pdf"),
  width = 8,
  height = 8
)

dba.plotPCA(dbObj,
            contrast=1,
            method=DBA_DESEQ2,
            attributes=DBA_CONDITION,
            label=DBA_ID
            )

dev.off()


# MA plot

pdf(
  file.path(output_dir, "Differential_enriched_MAplot.pdf"),
  width = 8,
  height = 8
)

dba.plotMA(dbObj,
           method=DBA_DESEQ2
           )

dev.off()

# Volcano plot

pdf(
  file.path(output_dir, "Differential_vol_plot.pdf"),
  width = 8,
  height = 8
)

dba.plotVolcano(dbObj,
                contrast = 1
                )

dev.off()

# Plot heatmap

pdf(
  file.path(output_dir, "Differential_enriched_heatmap.pdf"),
  width = 8,
  height = 8
)

hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
readscores <- dba.plotHeatmap(dbObj, correlations = FALSE,
                              scale="row", colScheme = hmap)

dev.off()

message(
  "Saved all differential peak plots "
)


### Relative Enrichment Results ###
# Extract results
res_deseq <- dba.report(res_deseq, method=DBA_DESEQ2, contrast = 1, th=1)

# Write GRanges to rds
saveRDS(res_deseq,
        file.path(output_dir, "res_deseq.rds"))

# Write results to a tsv file
out <- as.data.frame(res_deseq)
write.table(
  out,
  file = file.path(output_dir, "cKO_vs_WT_deseq2.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  out,
  file = file.path(output_dir, "cKO_vs_WT_deseq2.csv"),
  quote = FALSE,
  row.names = FALSE
)

# Create bed files for each keeping only significant peaks (p < 0.05)
cKO_enrich <- out %>% 
  filter(FDR < 0.05 & Fold > 0) %>% 
  dplyr::select(seqnames, start, end)
WT_enrich <- out %>% 
  filter(FDR < 0.05 & Fold < 0) %>% 
  dplyr::select(seqnames, start, end)

# Write to file
write.table(
  cKO_enrich,
  file = file.path(output_dir, "cKO_enriched.bed"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  WT_enrich,
  file = file.path(output_dir, "WT_enriched.bed"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(
  "Saved enriched peaks bed files, GRanges object, and relative enrichment files "
)

# NOTE: BED files cannot contain headers and so we have added the col.names=F argument
# to address that. Additionally, we took only the first three columns
# from the results (genomic coordinates) to adhere to a minimal BED file format.


## Debugging ##
# Read in samplesheet metadata
# output_dir <- "results/differential_enrichment"
# meta_csv <- read.csv("meta/diffbind_metadata.csv")
# samples <- names(meta_csv)

# Load in the existing DiffBind object saved on local
# dbObj <- readRDS("data/DiffBind/dbObj.rds")
#message("skipping dba count")
# Explore
#dbObj