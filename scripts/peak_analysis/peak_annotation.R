# Peak Annotation for ChIP using ChIPseeker

# example:
# Rscript peak_annotation.R \
# --narrowpeak_dir data/macs2/narrowPeak \
# --output_dir results/peak_annotation

# Load libraries
library(optparse)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)
library(ChIPpeakAnno)

#### INPUT ARGUMENT OPTIONS FOR COMMANDLINE #####

option_list <- list(
  make_option("--narrowpeak_dir", type = "character"),
  make_option("--output_dir", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))
output_dir <- opt$output_dir
narrowpeak_dir <- opt$narrowpeak_dir

# Argument Validations
if (is.null(narrowpeak_dir)) {
  stop("--narrowpeak_dir is required")
}

if (!dir.exists(narrowpeak_dir)) {
  stop(paste("Directory does not exist:", narrowpeak_dir))
}

if (is.null(output_dir)) {
  stop("--output_dir is required")
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# Load in peak files
sample_files <- list.files(
  path = narrowpeak_dir,
  pattern = "\\.narrowPeak$",
  full.names = TRUE
)

sample_names <- basename(sample_files)
sample_names <- sub("_peaks\\.narrowPeak$", "", sample_names)

if (length(sample_files) == 0) {
  stop(
    paste(
      "No narrowPeak files found in:",
      narrowpeak_dir
    )
  )
}

message("Found ", length(sample_files), " narrowPeak files")


granges_list <- list()

for (r in seq_along(sample_files)) {
  
  granges_list[[sample_names[r]]] <- ChIPpeakAnno::toGRanges(
    sample_files[r],
    format = "narrowPeak",
    header = FALSE
  )
}

# Load genome reference for annotation
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

### Peak Annotation ###

# Create a list of GRanges objects
samples_list <- granges_list

# Annotate each sample
peakAnnoList <- lapply(samples_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), annoDb="org.Mm.eg.db", verbose=FALSE)
# length(peakAnnoList) # Validate
# names(peakAnnoList) # Validate

# Loop to plot each annotated sample

for (sample_name in names(peakAnnoList)) {
  
  anno <- peakAnnoList[[sample_name]]
  
  # Save annotation table
  write.csv(
    as.data.frame(anno),
    file.path(
      output_dir,
      paste0(sample_name, "_annotation.csv")
    ),
    row.names = FALSE
  )
  
  # Annotation pie chart
  pdf(
    file.path(output_dir, paste0(sample_name, "_annotation_pie.pdf")),
    width = 8,
    height = 6
  )
  plotAnnoPie(anno)
  dev.off()
  
  # Annotation bar plot
  pdf(
    file.path(output_dir, paste0(sample_name, "_annotation_barplot.pdf")),
    width = 8,
    height = 6
  )
  plotAnnoBar(anno)
  dev.off()
  
  # UpSet plot
  pdf(
    file.path(output_dir, paste0(sample_name, "_annotation_upset.pdf")),
    width = 8,
    height = 6
  )
  upsetplot(anno)
  dev.off()
  
  # Distance to TSS
  pdf(
    file.path(output_dir, paste0(sample_name, "_distance_to_tss.pdf")),
    width = 8,
    height = 6
  )
  plotDistToTSS(anno)
  dev.off()
  
  message("Saved annotation plots for sample: ", sample_name)
}


## Visualization for multiple samples ##

# Plot barplot for all samples
pdf(
  file.path(output_dir, "all_samples_annotation_barplot.pdf"),
  width = 10,
  height = 8
)

plotAnnoBar(peakAnnoList)

dev.off()
message("Saved multi-sample annotation barplot")

# Plot distance to TSS for all samples
pdf(
  file.path(output_dir, "all_samples_distance_to_tss.pdf"),
  width = 10,
  height = 8
)

plotDistToTSS(peakAnnoList)

dev.off()
message("Saved multi-sample distance-to-TSS plot")

## Enrichment around TSS ##

# Get promoters
promoter <- getPromoters(TxDb = txdb, upstream = 2000, downstream = 2000)
promoter

# Create tag matrix (computationally expensive, deepTools_profile analog), optional
# for (sample_name in names(granges_list)) {
#   
#   tagMatrix <- getTagMatrix(
#     granges_list[[sample_name]],
#     windows = promoter
#   )
#   
#   pdf(
#     file.path(
#       output_dir,
#       paste0(sample_name, "_avg_profile.pdf")
#     ),
#     width = 8,
#     height = 6
#   )
#   
#   plotAvgProf(
#     tagMatrix,
#     xlim = c(-2000, 2000),
#     xlab = "Genomic Region (5'->3')",
#     ylab = "Read Count Frequency",
#     conf = 0.95,
#     resample = 1000
#   )
#   
#   dev.off()
  
  # heatmap profiles (requires a lot of memory!)
  # pdf(
  #   file.path(
  #     output_dir,
  #     paste0(sample_name, "_heatmap_profile.pdf")
  #   ),
  #   width = 8,
  #   height = 6
  # )
  # 
  # tagHeatmap(sample_name) # If you would like a heatmap profile for each sample
  # 
  # dev.off()
#}


## Plot Profiles for all samples (Requires larger RAM and CPU) ##

# Create tag matrix (computationally expensive, deepTools_profile analog), optional
# tagMatrixList <- lapply(
#   granges_list,
#   getTagMatrix,
#   windows = promoter
# )

# Plot profile plots, multiple lines in a single plot
#plotAvgProf(tagMatrixList, xlim=c(-2000, 2000))

# Plot profile plots, faceted by row
#plotAvgProf(tagMatrixList, xlim=c(-2000, 2000), conf=0.95,resample=500, facet="row")

# Plot multiple heatmaps (huge memory requirement), optional
#tagHeatmap(tagMatrixList)


## Debugging ##
# output_dir <- "results/peak_annotation" # Temporary to debug
# narrowpeak_dir <- "data/macs2/narrowPeak" # Temporary to debug

