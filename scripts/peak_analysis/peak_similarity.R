# Peak Similarity QC in R

# example executable:
# Rscript peak_similarity.R \
# --input_file data/multiBamSummary/multiBAMsummary_noInput.tab \
# --narrowpeak_dir data/macs2/narrowPeak \
# --output_dir results/peak_similarity

# This script requires a multiBamSummary.tab matrix

# Load Libraries
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(ggrepel)
library(ChIPpeakAnno)
library(UpSetR)
library(optparse)

#### INPUT ARGUMENT OPTIONS FOR COMMANDLINE #####

option_list <- list(
  make_option("--input_file", type = "character"),
  make_option("--output_dir", type = "character"),
  make_option("--narrowpeak_dir", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))
output_dir <- opt$output_dir
input_file <- opt$input_file
narrowpeak_dir <- opt$narrowpeak_dir

# Argument Validations
if (is.null(input_file)) {
  stop("--input_file is required")
}

if (is.null(output_dir)) {
  stop("--output_dir is required")
}

if (is.null(narrowpeak_dir)) {
  stop("--narrowpeak_dir is required")
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!dir.exists(narrowpeak_dir)) {
  stop(paste("Directory does not exist:", narrowpeak_dir))
}

# Read in Count Density

counts <- read.delim(input_file, sep = "\t")

# Requires multiBamSummary Matrix
# Here is an example script
#multiBamSummary bins \
#--bamfiles cKO_H3K27ac_ChIPseq_REP1.mLb.clN.sorted.bam cKO_H3K27ac_ChIPseq_REP2.mLb.clN.sorted.bam cKO_H3K27ac_ChIPseq_REP3.mLb.clN.sorted.bam \
#WT_H3K27ac_ChIPseq_REP1.mLb.clN.sorted.bam WT_H3K27ac_ChIPseq_REP2.mLb.clN.sorted.bam WT_H3K27ac_ChIPseq_REP3.mLb.clN.sorted.bam \
#--outFileName multiBamsummary_noInput.npz \
#--labels cKO_1 cKO_2 cKO_3 WT_1 WT_2 WT_3  \
#-p 6 \
#--outRawCounts multiBAMsummary_noInput.tab

# remove genomic coordinate info
plot_counts <- data.frame(counts[, 4:ncol(counts)])

# Change column names
colnames(plot_counts) <- colnames(plot_counts) %>% 
  str_replace( "X.", "") %>% 
  str_remove( "\\.$")

# Create meta
meta <- data.frame(row.names = colnames(plot_counts), 
                   genotype=colnames(plot_counts) %>% str_remove("\\_[0-9]"))

#### TRANSFORMING BAM COUNTS #####
# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(plot_counts, meta, design = ~genotype)

# Run vst and extract transformed counts
vst <- vst(dds)
vst_counts <- assay(vst)


#### PRINCIPAL COMPONENT ANALYSIS #####
# Compute principal components
pc <- prcomp(t(vst_counts))
plot_pca <- data.frame(pc$x, meta)
summary(pc) # will tell you how much variance is explained by each PC


# Plot with sample names used as data points
pca_plot <- ggplot(plot_pca, aes(PC1, PC2, color = genotype, label = rownames(plot_pca))) + 
  theme_bw() +
  geom_point() +
  geom_text_repel() +
  xlab("PC1") +
  ylab("PC2") +
  scale_x_continuous(expand = c(0.3,  0.3)) +
  theme(plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))

#print(pca_plot) # View plot

ggsave(
  file = file.path(output_dir, "PCA_plot.pdf" ),
  pca_plot,
  width = 8,
  height = 6
  )

message("Saved PCA plot")


#### INTER-SAMPLE CORRELATION #####

# Set annotation and colors
annotation <- meta
heat.colors <- brewer.pal(6, "YlOrRd")

# Plot ICA
pdf(
  file.path(output_dir, "Sample_correlation_heatmap.pdf"),
  width = 8,
  height = 8
)

pheatmap(
  cor(vst_counts),
  color = heat.colors,
  annotation = annotation
)

dev.off()

message(
  "Saved correlation heatmap: ",
  file.path(output_dir, "Sample_correlation_heatmap.pdf")
)


#### PEAK RANKING PLOT #####

# Get all narrowpeak file names and path
sample_files <- list.files(
  path = narrowpeak_dir,
  pattern = "\\.narrowPeak$",
  full.names = TRUE
)

# Create a vector of short names
vars <- basename(sample_files) %>%
  str_remove("_peaks\\.narrowPeak$")

# List of sample names
peak_list <- list()

# iterate
for(r in seq_along(sample_files)) {
  
  peaks <- read.delim(sample_files[r], header = FALSE)
  
  peak_list[[vars[r]]] <- data.frame(
    peak_enrichment = peaks$V7,
    peak_rank = rank(dplyr::desc(peaks$V7))
  ) %>%
    arrange(peak_rank)
}


#  WT only
wt <- bind_rows(
  peak_list[c(
    "WT_H3K27ac_ChIPseq_REP1",
    "WT_H3K27ac_ChIPseq_REP2",
    "WT_H3K27ac_ChIPseq_REP3"
  )],
  .id = "reps"
)

wt_rank_plot <- ggplot(
  wt,
  aes(peak_rank, peak_enrichment, color = reps)
) +
  geom_line() +
  ggtitle("WT samples") +
  theme_bw() +
  geom_hline(yintercept = 4, linetype = "dashed", color = "red") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  xlab("Peak rank") +
  ylab("Peak enrichment")

print(wt_rank_plot)

ggsave(
  filename = file.path(output_dir, "WT_peak_rank_plot.pdf"),
  plot = wt_rank_plot,
  width = 8,
  height = 6
)

message(
  "Saved WT ranked peaks plot: ",
  file.path(output_dir, "WT_peak_rank_plot.pdf")
)

# cKO only
cko <- bind_rows(
  peak_list[c(
    "cKO_H3K27ac_ChIPseq_REP1",
    "cKO_H3K27ac_ChIPseq_REP2",
    "cKO_H3K27ac_ChIPseq_REP3"
  )],
  .id = "reps"
)

cko_rank_plot <- ggplot(
  cko,
  aes(peak_rank, peak_enrichment, color = reps)
) +
  geom_line() +
  ggtitle("cKO samples") +
  theme_bw() +
  geom_hline(yintercept = 4, linetype = "dashed", color = "red") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  xlab("Peak rank") +
  ylab("Peak enrichment")

print(cko_rank_plot)

ggsave(
  filename = file.path(output_dir, "cko_peak_rank_plot.pdf"),
  plot = cko_rank_plot,
  width = 8,
  height = 6
)

message(
  "Saved cko ranked peaks plot: ",
  file.path(output_dir, "cko_peak_rank_plot.pdf")
)


# Combine data into one dataframe
allreps <- rbind(wt, cko)

# Add a column for genotype
allreps$genotype <- str_split_fixed(allreps$reps, "_", 4)[,1]

# Plot histogram
chip_hist_plot <- ggplot(allreps, aes(peak_enrichment, fill = genotype)) +
  geom_histogram(aes(peak_enrichment), binwidth = 3) +
  ggtitle("Histogram of peak enrichment values") +
  theme_bw() +
  geom_vline(xintercept=median(allreps$peak_enrichment), linetype='dashed') +
  geom_vline(xintercept=4, linetype='dashed', color="red") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Peak Ranked Enrichment") 
  print(chip_hist_plot)

ggsave(
  filename = file.path(output_dir, "peak_enrichment_histogram.pdf"),
  plot = chip_hist_plot,
  width = 8,
  height = 6
)

message(
  "Saved histogram enrichment plot: ",
  file.path(output_dir, "peak_enrichment_histogram.pdf")
)


#### PEAK OVERLAPPING CONSISTENCY CHECK #####

# Load libraries
library(GenomicRanges)
library(IRanges)

# Loop across list
granges_list <- list()

for(r in seq_along(sample_files)) {
  
  granges_list[[vars[r]]] <- ChIPpeakAnno::toGRanges(
    sample_files[r],
    format = "narrowPeak",
    header = FALSE
  )
}

# print(names(granges_list)) #validation check

# Find overlapping peaks and merge
# connectedPeaks to "merge" will add 1 to the overlapping counts

# alias variables are required for ChIPpeakAnno to correctly assign in the plot
WT_REP1 <- granges_list[["WT_H3K27ac_ChIPseq_REP1"]]
WT_REP2 <- granges_list[["WT_H3K27ac_ChIPseq_REP2"]]
WT_REP3 <- granges_list[["WT_H3K27ac_ChIPseq_REP3"]]

# NOTE:
# findOverlapsOfPeaks() may emit Seqinfo warnings when individual
# narrowPeak files contain peaks on different random/unplaced contigs.
# This does not necessarily indicate different genome builds.
# for now I have suppressed this warning.

suppressWarnings({ # Be careful to check objects based on the same reference genome
  olaps_wt <- findOverlapsOfPeaks(
    WT_REP1,
    WT_REP2,
    WT_REP3,
    connectedPeaks = "merge")
})

# class(olaps_wt) #validation check
# names(olaps_wt) #validation check
# colnames(olaps_wt$venn_cnt) #validation check and fighting the package

# Venn Diagram of Peak Overlaps / Concordance
pdf(
  file.path(output_dir, "WT_peak_overlaps_venndiagram.pdf"),
  width = 8,
  height = 8
)

# plot Venn Diagram

# NOTE:
# makeVennDiagram() performs a HyperG overlap calculation.
# No totalTest universe is supplied because this figure is used
# for replicate concordance QC rather than overlap significance testing.
# The resulting warning can be ignored.

suppressWarnings({
  venstats <- makeVennDiagram(
    olaps_wt,
    connectedPeaks = "merge",
    fill = c("#CC79A7", "#56B4E9", "#F0E442"),
    col = c("#D55E00", "#0072B2", "#E69F00"),
    cat.col = c("#D55E00", "#0072B2", "#E69F00")
  )
})

dev.off()
message("Saved WT overlap Venn diagram")

# Prepare data for UpSetR
pdf(
  file.path(output_dir, "WT_peak_overlap_upsetplot.pdf"),
  width = 9,
  height = 7,
  onefile = FALSE # prevents blank first page
)

set_counts <- olaps_wt$venn_cnt[, colnames(olaps_wt$venn_cnt)] %>%
  as.data.frame() %>%
  mutate(group_number = row_number()) %>%
  pivot_longer(!Counts & !group_number, names_to = 'sample', values_to = 'member') %>%
  filter(member > 0) %>%
  group_by(Counts, group_number) %>%
  summarize(group = paste(sample, collapse = "&"), .groups = "drop")

# Set required variables
set_counts_upset <- set_counts$Counts
names(set_counts_upset) <- set_counts$group

# Plot the UpSet plot
upset(fromExpression(set_counts_upset), order.by = "freq", text.scale = 1.3)

dev.off()

message(
  "Saved UpSet plot: ",
  file.path(output_dir, "WT_peak_overlap_upsetplot.pdf")
)


## Debugging ##
# counts <- read.delim("data/multiBamSummary/multiBAMsummary_noinput.tab", sep = "\t") # Temporary to debug
# output_dir <- "results/peak_similarity" # Temporary to debug
# narrowpeak_dir <- "data/macs2/narrowPeak" # Temporary to debug