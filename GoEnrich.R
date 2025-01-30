setwd("/Users/zhouxin/Downloads/RNAseq_course")

# Load necessary packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(cowplot)
library(patchwork)

# Function: Perform GO enrichment analysis and save results
run_go_enrichment <- function(de_genes, all_genes, comparison_name, gene_list) {
  if (length(de_genes) == 0) {
    warning(paste0("No significant genes found for ", comparison_name))
    return(NULL)
  }
  
  # Perform GO enrichment analysis
  go_enrichment <- enrichGO(
    gene = de_genes,
    universe = all_genes,
    OrgDb = org.Hs.eg.db,
    ont = "BP",  # Biological process
    keyType = "ENSEMBL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  # Save enrichment results as a CSV file
  write.csv(go_enrichment@result, file = paste0("GO_Enrichment_", comparison_name, ".csv"))
  
  # Generate a dot plot
  dotplot(go_enrichment, showCategory = 30, font.size = 14) +  # Increase font size
    ggtitle(paste0("GO Enrichment: ", comparison_name))
  ggsave(filename = paste0("GO_Dotplot_", comparison_name, ".jpeg"), width = 16, height = 12, dpi = 300)
  
  # Generate a bar plot
  barplot(go_enrichment, showCategory = 20, title = paste0("GO Enrichment: ", comparison_name))
  ggsave(filename = paste0("GO_Barplot_", comparison_name, ".jpeg"), width = 16, height = 12, dpi = 300)
  
  # Convert gene IDs to readable format
  edox <- setReadable(go_enrichment, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
  
  # Gene-Concept Network plot (showing primary relationships only)
  cnet_p1 <- cnetplot(edox, foldChange = gene_list)
  ggsave(filename = paste0("GO_Cnetplot_", comparison_name, ".jpeg"), plot = cnet_p1, width = 12, height = 12, dpi = 300)
  
  # Heatplot (showing primary relationships only)
  heat_p1 <- heatplot(edox, showCategory = 10, foldChange = gene_list)  # Increase displayed categories
  ggsave(filename = paste0("GO_Heatplot_", comparison_name, ".jpeg"), plot = heat_p1, width = 16, height = 10, dpi = 300)
  
  # Treeplot (showing primary relationships only)
  edox2 <- pairwise_termsim(edox)
  tree_p1 <- treeplot(edox2)
  ggsave(filename = paste0("GO_Treeplot_", comparison_name, ".jpeg"), plot = tree_p1, width = 12, height = 12, dpi = 300)
  
  print(paste0("Top 10 GO terms for ", comparison_name, ":"))
  print(head(go_enrichment@result, 10))
  
  return(go_enrichment)
}

# Define the background gene set
background_genes <- rownames(dds)

# TNBC vs Normal
res_tnbc_vs_normal <- results(dds, contrast = c("condition", "TNBC", "Normal"))
de_genes_tnbc <- rownames(subset(res_tnbc_vs_normal, padj < 0.05))
gene_list_tnbc <- res_tnbc_vs_normal$log2FoldChange
names(gene_list_tnbc) <- rownames(res_tnbc_vs_normal)
go_tnbc <- run_go_enrichment(de_genes_tnbc, background_genes, "TNBC_vs_Normal", gene_list_tnbc)

# HER2 vs Normal
res_her2_vs_normal <- results(dds, contrast = c("condition", "HER2", "Normal"))
de_genes_her2 <- rownames(subset(res_her2_vs_normal, padj < 0.05))
gene_list_her2 <- res_her2_vs_normal$log2FoldChange
names(gene_list_her2) <- rownames(res_her2_vs_normal)
go_her2 <- run_go_enrichment(de_genes_her2, background_genes, "HER2_vs_Normal", gene_list_her2)

# NonTNBC vs Normal
res_nontnbc_vs_normal <- results(dds, contrast = c("condition", "NonTNBC", "Normal"))
de_genes_nontnbc <- rownames(subset(res_nontnbc_vs_normal, padj < 0.05))
gene_list_nontnbc <- res_nontnbc_vs_normal$log2FoldChange
names(gene_list_nontnbc) <- rownames(res_nontnbc_vs_normal)
go_nontnbc <- run_go_enrichment(de_genes_nontnbc, background_genes, "NonTNBC_vs_Normal", gene_list_nontnbc)
