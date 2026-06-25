# Differential peaks analysis

# example:
# Rscript differential_peaks_analysis.R \
# --input_dir results/differential_enrichment \
# --output_dir results/differential_peak_analysis

# Load libraries 
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(optparse)

## DEBUG ###
# Read in samplesheet metadata
output_dir <- "results/differential_peak_analysis"
input_dir <- "results/differential_enrichment"


# Load result from DiffBind analysis if not present in your environment
res_deseq <- readRDS("results/differential_enrichment/res_deseq.rds")

# Set the annotation database if not present in your environment
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


# Add annotations to our DiffBind results
annot_res_all <- annotatePeak(res_deseq,
                              tssRegion = c(-3000, 3000),
                              TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                              annoDb = "org.Mm.eg.db")

# Barplot
pdf(
  file.path(output_dir, "de_annotated_barplot.pdf"),
  width = 6,
  height = 4
)

plotAnnoBar(annot_res_all)

dev.off()

# TSS distance plot
pdf(
  file.path(output_dir, "de_TSSdistance_plot.pdf"),
  width = 6,
  height = 4
)

plotDistToTSS(annot_res_all)

dev.off()

## FUNCTIONAL ANALYSIS VIA CLUSTERPROFILER ##

# Create a dataframe from anno results

annot_res_all_df <- as.data.frame(annot_res_all)

# Save a table for all de peaks

write.table(
  annot_res_all_df,
  file = file.path(output_dir, "annotated_de_peaks_cko_wt.xls"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Get background gene set of Entrez IDs
background_set <- as.character(annot_res_all_df$geneId)

# Prepare gene set query for up-regulated genes
sigUp <- dplyr::filter(annot_res_all_df, FDR < 0.05, Fold > 0)
sigUp_genes <- as.character(sigUp$geneId)

# Run over-representation analysis
go_ORA_Up <- enrichGO(gene = sigUp_genes,
                      universe = background_set,
                      keyType = "ENTREZID",
                      OrgDb = org.Mm.eg.db,
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)

# Save data frame to file
go_ORA_Up_df <- data.frame(go_ORA_Up)
write.csv(go_ORA_Up_df,
          file = file.path(output_dir, "GO_ORA_clusterProfiler_cKO_vs_WT_Upregulated.csv")
          )
         
View(go_ORA_Up_df)


# Dotplot
pdf(
  file.path(output_dir, "de_ORA_dotplot.pdf"),
  width = 6,
  height = 8
)

dotplot(go_ORA_Up)

dev.off()

# Before generating the enrichment plot, the similarity between 
# terms must be computed using the pairwise_termsim() function.


# Enrich plot
go_ORA_Up <- enrichplot::pairwise_termsim(go_ORA_Up)

pdf(
  file.path(output_dir, "de_ORA_enrichment_plot.pdf"),
  width = 6,
  height = 8
)

emapplot(go_ORA_Up)

dev.off()

# Lastly use GREAT comparative analysis for distal enriched regions
# Also check MEME, DAVID, clusterProfiler, and ReviGO for more!