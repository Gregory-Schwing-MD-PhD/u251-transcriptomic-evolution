#!/usr/bin/env Rscript

# ==============================================================================
# U251 EVOLUTION & KITCHEN SINK VISUALIZATION SUITE
# ==============================================================================

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(dplyr)
    library(RColorBrewer)
    library(ape)      # For Phylogenetic trees
    library(ggrepel)  # For non-overlapping labels
})

args <- commandArgs(trailingOnly = TRUE)
results_dir   <- args[1]
output_prefix <- args[2]
counts_file   <- args[3]
meta_file     <- args[4]
contrast_file <- args[5]

dir.create(dirname(output_prefix), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. DATA LOADING & SETUP
# ==============================================================================
cat("LOG: Loading Data...\n")

# Load Metadata
meta <- read.csv(meta_file, row.names=1)
# Enforce Biological Order: Culture -> Primary -> Recurrent
meta$Classification <- factor(meta$Classification, levels = c("Culture_U2", "Primary_U2", "Recurrent_U2"))

# Load Counts
counts <- read.table(counts_file, header=TRUE, row.names=1, check.names=FALSE)
counts <- counts[, rownames(meta)] # Sync samples

# Create DESeq Object
dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = meta, design = ~ Classification)

# ==============================================================================
# 2. EVOLUTIONARY ANALYSIS (LRT & TREES)
# ==============================================================================
cat("LOG: Running Likelihood Ratio Test (LRT) for Evolution...\n")

# Run LRT (changes across ANY timepoint)
dds_lrt <- DESeq(dds, test="LRT", reduced=~1)
res_lrt <- results(dds_lrt)
vsd <- vst(dds_lrt, blind=FALSE)
mat_vst <- assay(vsd)

# A. ROOTED PHYLOGENETIC TREE
cat("LOG: Generating Rooted Phylogenetic Tree...\n")
sig_genes_lrt <- rownames(res_lrt[which(res_lrt$padj < 0.05), ])
mat_sig <- mat_vst[sig_genes_lrt, ]
dists <- dist(t(mat_sig))
hc <- hclust(dists, method = "ward.D2")
phylo_tree <- as.phylo(hc)

# Try to root at Culture
tryCatch({
    rooted_tree <- root(phylo_tree, outgroup = "C2B", resolve.root = TRUE)
    pdf(paste0(output_prefix, "_Phylogenetic_Tree.pdf"), width=8, height=6)
    plot(rooted_tree, main="Phylogenetic Reconstruction of U251 Evolution", type="phylogram", edge.width=2, label.offset=0.5, cex=1.2)
    add.scale.bar()
    tiplabels(pch=21, bg=c("blue", rep("orange", 3), rep("red", 3)), cex=2)
    legend("topleft", legend=c("Culture", "Primary", "Recurrent"), pch=21, pt.bg=c("blue", "orange", "red"))
    dev.off()
}, error = function(e) { cat("WARN: Could not root tree at 'C2B'. Check sample names.\n") })

# B. EVOLUTIONARY TRAJECTORY (PCA)
cat("LOG: Generating PCA Trajectory...\n")
pca <- prcomp(t(mat_sig))
pcaData <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Sample = rownames(meta), Class = meta$Classification)
centroids <- aggregate(cbind(PC1, PC2) ~ Class, data=pcaData, mean)

p_pca <- ggplot(pcaData, aes(x=PC1, y=PC2, color=Class)) +
  geom_point(size=5) +
  geom_path(data=centroids, aes(group=1), arrow=arrow(length=unit(0.3, "cm"), type="closed"), color="black", alpha=0.5, linetype="dashed") +
  geom_text_repel(aes(label=Sample)) +
  scale_color_manual(values=c("Culture_U2"="blue", "Primary_U2"="orange", "Recurrent_U2"="red")) +
  theme_bw() + ggtitle("Evolutionary Trajectory (PCA)")
ggsave(paste0(output_prefix, "_PCA_Trajectory.pdf"), p_pca, width=8, height=6)

# ==============================================================================
# 3. KITCHEN SINK (STANDARD PLOTS)
# ==============================================================================
cat("LOG: Generating Standard QC Plots...\n")

# C. SAMPLE DISTANCE HEATMAP
sampleDists <- dist(t(mat_vst))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf(paste0(output_prefix, "_SampleDistances.pdf"), width=8, height=8)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, main="Sample-to-Sample Distances")
dev.off()

# D. TOP 50 GENES HEATMAP
top_genes <- head(order(res_lrt$padj), 50)
mat_top <- mat_vst[top_genes, ] - rowMeans(mat_vst[top_genes, ])
pdf(paste0(output_prefix, "_Heatmap_Top50_Evolution.pdf"), width=8, height=10)
pheatmap(mat_top, annotation_col=as.data.frame(colData(dds_lrt)[,c("Classification","Environment")]), show_rownames=TRUE, cluster_cols=FALSE, main="Top 50 Evolution Genes")
dev.off()

# ==============================================================================
# 4. AUTOMATED VOLCANO PLOTS (FROM CONTRASTS.CSV)
# ==============================================================================
cat("LOG: Generating Pairwise Volcano Plots...\n")

# Re-run DESeq with Wald test for pairwise comparisons
dds_wald <- DESeq(dds, test="Wald")
contrasts <- read.csv(contrast_file)

for(i in 1:nrow(contrasts)) {
    ref_group <- contrasts$reference[i]
    target_group <- contrasts$target[i]
    comp_name <- contrasts$id[i]
    
    cat(sprintf("   -> Processing Contrast: %s vs %s (%s)\n", ref_group, target_group, comp_name))
    
    res <- results(dds_wald, contrast=c("Classification", target_group, ref_group))
    
    # Prepare data for plotting
    res_df <- as.data.frame(res)
    res_df$sig <- "NS"
    res_df$sig[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "UP"
    res_df$sig[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "DOWN"
    
    p_vol <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=sig)) +
        geom_point(alpha=0.5) +
        scale_color_manual(values=c("DOWN"="blue", "NS"="grey", "UP"="red")) +
        theme_bw() +
        geom_vline(xintercept=c(-1, 1), linetype="dashed") +
        geom_hline(yintercept=-log10(0.05), linetype="dashed") +
        ggtitle(paste0("Volcano: ", target_group, " vs ", ref_group))
        
    ggsave(paste0(output_prefix, "_Volcano_", comp_name, ".pdf"), p_vol, width=6, height=5)
}

cat("LOG: Analysis Complete.\n")
