# Peak Similarity QC in R

# Load Libraries
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(ggrepel)
library(ChIPpeakAnno)
library(UpSetR)

# Read Count Density
# Read in data
counts <- read.delim("data/multiBamSummary/multiBAMsummary_noinput.tab", sep="\t")

# Command Line script to generate matrix
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

# Transforming the BAM Counts
# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(plot_counts, meta, design = ~genotype)

# Run vst and extract transformed counts
vst <- vst(dds)
vst_counts <- assay(vst)

# Compute principal components
pc <- prcomp(t(vst_counts))
plot_pca <- data.frame(pc$x, meta)
summary(pc) # will tell you how much variance is explained by each PC

# Plot with sample names used as data points
ggplot(plot_pca, aes(PC1, PC2, color = genotype, label = rownames(plot_pca)), size = 3 ) + 
  theme_bw() +
  geom_point() +
  geom_text_repel() +
  xlab('PC1 (33% of variance)') +
  ylab('PC2 (22% of variance)') +
  scale_x_continuous(expand = c(0.3,  0.3)) +
  theme(plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))


# Inter-sample correlation

# Set annotation and colors
annotation <- meta
heat.colors <- brewer.pal(6, "YlOrRd")

# Plot ICA
pheatmap(cor(vst_counts), color=heat.colors, annotation=annotation)


# Peak Ranking Plot
# Get all narrowpeak file names and path
sample_files <- list.files(path = "./data/macs2/narrowPeak/", full.names = T)

# Create a vector of short names
vars <- str_remove( sample_files, "./data/macs2/narrowPeak/") %>% 
  str_remove("_peaks.narrowPeak")

# Loop through to create a dataframe for each sample with columns required
for(r in 1:length(sample_files)){
  peaks <- read.delim(sample_files[r], header = FALSE)
  df <- data.frame(peak_enrichment = peaks$V7, peak_rank = rank(dplyr::desc(peaks$V7))) %>% 
    dplyr::arrange(peak_rank) 
  assign(vars[r], df)
}

#  WT only
wt <- bind_rows("WT_H3K27ac_ChIPseq_REP1" = WT_H3K27ac_ChIPseq_REP1, 
                "WT_H3K27ac_ChIPseq_REP2" = WT_H3K27ac_ChIPseq_REP2,
                "WT_H3K27ac_ChIPseq_REP3" = WT_H3K27ac_ChIPseq_REP3,
                .id = "reps")

ggplot(wt, aes(peak_rank, peak_enrichment, color = reps)) + 
  geom_line() +
  ggtitle("WT samples") +
  theme_bw() +
  geom_hline(yintercept=4, linetype='dashed', color="red") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Peak rank") + ylab("Peak enrichment")


# cKO only
cko <- bind_rows("cKO_H3K27ac_ChIPseq_REP1" = cKO_H3K27ac_ChIPseq_REP1, 
                 "cKO_H3K27ac_ChIPseq_REP2" = cKO_H3K27ac_ChIPseq_REP2,
                 "cKO_H3K27ac_ChIPseq_REP3" = cKO_H3K27ac_ChIPseq_REP3,
                 .id = "reps")

ggplot(cko, aes(peak_rank, peak_enrichment, color = reps)) + 
  geom_line() +
  ggtitle("cKO samples") +
  theme_bw() +
  geom_hline(yintercept=4, linetype='dashed', color="red") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Peak rank") + ylab("Peak enrichment")

# Combine data into one dataframe
allreps <- rbind(wt, cko)

# Add a column for genotype
allreps$genotype <- str_split_fixed(allreps$reps, "_", 4)[,1]

# Plot histogram
ggplot(allreps, aes(peak_enrichment, fill = genotype)) + 
  geom_histogram(aes(peak_enrichment)) +
  ggtitle("Histogram of peak enrichment values") +
  theme_bw() +
  geom_vline(xintercept=median(df$peak_enrichment), linetype='dashed') +
  geom_vline(xintercept=4, linetype='dashed', color="red") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Peak enrichment")

#### PEAK OVERLAPPING CONSISTENCY CHECK #####

library(GenomicRanges)
library(IRanges)

for(r in 1:length(sample_files)){
  obj <- ChIPpeakAnno::toGRanges(sample_files[r], format="narrowPeak", header=FALSE)  
  assign(sample_names[r], obj)
}

# Find overlapping peaks and merge
# connectedPeaks to "merge" will add 1 to the overlapping counts

olaps_wt <- findOverlapsOfPeaks(WT_H3K27ac_ChIPseq_REP1,
                                WT_H3K27ac_ChIPseq_REP2,
                                WT_H3K27ac_ChIPseq_REP3, connectedPeaks = "merge")

# Venn Diagram of Peak Overlaps / Concordance
venstats <- makeVennDiagram(olaps_wt, connectedPeaks = "merge",
                            fill=c("#CC79A7", "#56B4E9", "#F0E442"), # circle fill color
                            col=c("#D55E00", "#0072B2", "#E69F00"), #circle border color
                            cat.col=c("#D55E00", "#0072B2", "#E69F00")) # category name color

# Prepare data for UpSetR
set_counts <- olaps_wt$venn_cnt[, colnames(olaps_wt$venn_cnt)] %>% 
  as.data.frame() %>% 
  mutate(group_number = row_number()) %>%
  pivot_longer(!Counts & !group_number, names_to = 'sample', values_to = 'member') %>%
  filter(member > 0) %>%
  group_by(Counts, group_number) %>% 
  summarize(group = paste(sample, collapse = '&'))

# Set required variables 
set_counts_upset <- set_counts$Counts
names(set_counts_upset) <- set_counts$group

# Plot the UpSet plot
upset(fromExpression(set_counts_upset), order.by = "freq", text.scale = 1.5)

