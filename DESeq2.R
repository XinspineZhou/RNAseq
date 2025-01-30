# Check and install required R packages
required_packages <- c("DESeq2", "ggplot2", "pheatmap", "BiocManager")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (pkg == "BiocManager") {
      install.packages("BiocManager")
    } else {
      BiocManager::install(pkg)
    }
  }
}
BiocManager::install("DESeq2")
if (!require("apeglm", quietly = TRUE)) {
  BiocManager::install("apeglm")
}

# Load installed packages
library(DESeq2)
library(ggplot2)
library(pheatmap)

# Define sample information table
sample_info <- data.frame(
  sample = c("TNBC1", "TNBC2", "TNBC3",
             "Normal1", "Normal2", "Normal3",
             "NonTNBC1", "NonTNBC2", "NonTNBC3",
             "HER21", "HER22", "HER23"),
  condition = c("TNBC", "TNBC", "TNBC",
                "Normal", "Normal", "Normal",
                "NonTNBC", "NonTNBC", "NonTNBC",
                "HER2", "HER2", "HER2"),
  replicate = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3) # Biological replicates
)

rownames(sample_info) <- sample_info$sample
print(sample_info)

# Set working directory (update with your actual path)
setwd("/Users/zhouxin/Downloads/RNAseq_course")

# List all _counts.txt files
file_list <- list.files(pattern = "_counts\\.txt$")
print(file_list)  # Check listed files

# Function: Read each count file and extract GeneID & count column
read_count_file <- function(file) {
  data <- read.table(file, header = TRUE, row.names = 1, sep = "\t")
  sample_name <- gsub("_counts\\.txt$", "", file)
  counts <- data[, ncol(data), drop = FALSE]  # Extract last column (count values)
  colnames(counts) <- sample_name
  return(counts)
}

# Read all files and combine them into a count matrix
all_counts <- lapply(file_list, read_count_file)
count_matrix <- do.call(cbind, all_counts)
head(count_matrix)

# Save the count matrix
write.csv(count_matrix, "combined_counts_matrix.csv")

# Ensure columns match sample info row names
if (!all(colnames(count_matrix) %in% rownames(sample_info))) {
  stop("Mismatch between count matrix columns and sample info rows!")
}

# Reorder count matrix columns to match sample info
count_matrix <- count_matrix[, rownames(sample_info)]
if (!all(colnames(count_matrix) == rownames(sample_info))) {
  stop("Sample names in count matrix and sample info are not in the same order!")
}

# Step 2: Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition
)

# Filter low-expression genes (genes with >=10 counts in at least 3 samples)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# Run DESeq2 analysis
dds <- DESeq(dds)

# Step 3: Extract differential expression results
res <- results(dds)
res <- lfcShrink(dds, coef = 2, type = "apeglm")

# Save differential expression results
write.csv(as.data.frame(res), "differential_expression_results.csv")

# Step 4: Visualization

# Volcano plot
res$significant <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "yes", "no")
volcano_plot <- ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  theme_minimal() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value")
ggsave("Volcano_plot.jpeg", plot = volcano_plot, width = 6, height = 4)

# PCA plot
vsd <- vst(dds, blind = TRUE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()
ggsave("PCA_plot_condition.jpeg", plot = pca_plot, width = 6, height = 4)

# Heatmap (Top 20 most variable genes)
top_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
pheatmap(assay(vsd)[top_genes, ], cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, annotation_col = as.data.frame(colData(dds)["condition"]),
         filename = "Heatmap_top20_genes.jpeg")

# Step 5: Save normalized counts
normalized_counts <- as.data.frame(assay(vsd))
write.csv(normalized_counts, "normalized_counts_matrix.csv")

# Step 6: Sample-to-sample distance heatmap
sample_dist <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dist)
pheatmap(sample_dist_matrix, clustering_distance_rows = sample_dist, clustering_distance_cols = sample_dist,
         annotation_col = as.data.frame(colData(dds)["condition"]), main = "Sample-to-Sample Distance Heatmap",
         filename = "Sample_to_Sample_Heatmap.jpeg")

# Step 7: Differential expression analysis for specific comparisons

# Load annotation database
if (!require("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)

# Define genes of interest
genes_of_interest <- c("ENSG00000113140", "ENSG00000204628", "ENSG00000127022")  # SPARC, RACK1, CANX

# Function: Perform DE analysis and visualization
run_differential_expression_analysis <- function(comparison_name, condition1, condition2, dds, genes_of_interest) {
  res <- results(dds, contrast = c("condition", condition1, condition2))
  significant_genes <- subset(res, padj < 0.05)
  
  num_upregulated <- sum(significant_genes$log2FoldChange > 0)
  num_downregulated <- sum(significant_genes$log2FoldChange < 0)
  
  message(comparison_name, " Total DE genes: ", nrow(significant_genes))
  message(comparison_name, " Upregulated: ", num_upregulated)
  message(comparison_name, " Downregulated: ", num_downregulated)
  
  write.csv(as.data.frame(significant_genes), file = paste0(comparison_name, "_significant_genes.csv"))

  res$gene_highlight <- ifelse(rownames(res) %in% genes_of_interest, "interest",
                               ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "yes", "no"))
  
  volcano_plot <- ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = gene_highlight), alpha = 0.6, size = 2) +
    scale_color_manual(values = c("interest" = "green", "yes" = "red", "no" = "grey")) +
    theme_minimal() +
    ggtitle(paste0("Volcano Plot: ", comparison_name))
  
  ggsave(paste0("Volcano_plot_", comparison_name, ".jpeg"), plot = volcano_plot, width = 8, height = 6)
  
  return(significant_genes)
}

# Run DE analysis for different comparisons
run_differential_expression_analysis("NonTNBC_vs_Normal", "NonTNBC", "Normal", dds, genes_of_interest)
run_differential_expression_analysis("HER2_vs_Normal", "HER2", "Normal", dds, genes_of_interest)
run_differential_expression_analysis("TNBC_vs_Normal", "TNBC", "Normal", dds, genes_of_interest)
