# Peak Annotation of Differentially Enriched Peaks

# Load libraries 
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)

# Load result from DiffBind analysis if not present in your environment
res_deseq <- readRDS("all_res_deseq2.rds")

# Set the annotation database if not present in your environment
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


# Add annotations to our DiffBind results
annot_res_all <- annotatePeak(res_deseq,
                              tssRegion = c(-3000, 3000),
                              TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                              annoDb = "org.Mm.eg.db")

# Barplot
plotAnnoBar(annot_res_all)

# TSS distance plot
plotDistToTSS(annot_res_all)


## FUNCTIONAL ANALYSIS VIA CLUSTERPROFILER ##

# Create a dataframe from anno results
annot_res_all_df <- as.data.frame(annot_res_all)

# Save a table for all de peaks
write.table(annot_res_all_df, file="results/annotated_de_peaks_cko_wt.xls", sep="\t", quote=F, row.names=F)

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
write.csv(go_ORA_Up_df, "results/GO_ORA_clusterProfiler_cKO_vs_WT_Upregulated.csv")
View(go_ORA_Up_df)


# Dotplot
dotplot(go_ORA_Up)

# Before generating the enrichment plot, the similarity between 
# terms must be computed using the pairwise_termsim() function.


# Enrich plot
go_ORA_Up <- enrichplot::pairwise_termsim(go_ORA_Up)
emapplot(go_ORA_Up)