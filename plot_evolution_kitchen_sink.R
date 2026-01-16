#!/usr/bin/env Rscript

# ==============================================================================
# U251 EVOLUTION & KITCHEN SINK VISUALIZATION SUITE (MEGA EDITION v2.1)
# Features:
# 1. Phylogenetic Trees & Trajectories (LRT)
# 2. Global GSVA (ComplexHeatmap)
# 3. Pairwise EnhancedVolcano + GSEA
# UPDATES:
# - Swapped org.Hs.eg.db for EnsDb.Hsapiens.v86 for better ID coverage.
# ==============================================================================

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(dplyr)
    library(RColorBrewer)
    library(ape)
    library(ggrepel)
    # library(org.Hs.eg.db) # <-- REMOVED
    library(EnsDb.Hsapiens.v86) # <-- ADDED
    library(clusterProfiler)
    library(enrichplot)
    library(GSVA)
    library(GSEABase)
    library(ComplexHeatmap)
    library(circlize)
    library(EnhancedVolcano)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) stop("Usage: script.R <results_dir> <out_prefix> <counts> <meta> <contrasts> <gmt>")

results_dir   <- args[1]
output_prefix <- args[2]
counts_file   <- args[3]
meta_file     <- args[4]
contrast_file <- args[5]
gmt_file      <- args[6]

dir.create(dirname(output_prefix), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. DATA LOADING & ID MAPPING
# ==============================================================================
cat("LOG: Loading Data...\n")
meta <- read.csv(meta_file, row.names=1)
# Ensure factors match your experimental design
meta$Classification <- factor(meta$Classification, levels = c("Culture_U2", "Primary_U2", "Recurrent_U2"))

counts <- read.table(counts_file, header=TRUE, row.names=1, check.names=FALSE)
counts <- counts[, rownames(meta)]

cat("LOG: Mapping IDs to Symbols (Ensembl)...\n")
clean_ids <- sub("\\..*", "", rownames(counts))

# Map using EnsDb (Native Ensembl Support)
syms <- suppressWarnings(mapIds(EnsDb.Hsapiens.v86, 
                                keys=clean_ids, 
                                column="SYMBOL", 
                                keytype="GENEID", 
                                multiVals="first"))

# Fallback for NAs
syms[is.na(syms)] <- clean_ids[is.na(syms)]
rownames(counts) <- make.unique(as.character(syms))

dds <- DESeqDataSetFromMatrix(round(counts), meta, ~ Classification)

# ==============================================================================
# 2. EVOLUTIONARY ANALYSIS (LRT)
# ==============================================================================
cat("LOG: Running Evolution LRT...\n")
dds_lrt <- DESeq(dds, test="LRT", reduced=~1)
res_lrt <- results(dds_lrt)
vsd <- vst(dds_lrt, blind=FALSE)
mat_vst <- assay(vsd)

# --- A. Phylogenetic Tree ---
cat("LOG: generating Tree...\n")
sig_genes <- rownames(res_lrt[which(res_lrt$padj < 0.05), ])
mat_sig <- mat_vst[sig_genes, ]
dists <- dist(t(mat_sig))
phylo_tree <- as.phylo(hclust(dists, method="ward.D2"))

tryCatch({
    # Attempt to root at C2B (Culture control) if it exists in your sample names
    rooted_tree <- root(phylo_tree, outgroup="C2B", resolve.root=TRUE)
    pdf(paste0(output_prefix, "_Phylogenetic_Tree.pdf"), width=8, height=6)
    plot(rooted_tree, main="Phylogenetic Reconstruction", type="phylogram", edge.width=2, cex=1.2)
    tiplabels(pch=21, bg=c("blue", rep("orange", 3), rep("red", 3)), cex=2)
    dev.off()
}, error = function(e) {
    cat("WARN: Could not root tree at C2B (Sample not found or other error).\n")
    # Plot unrooted if rooting fails
    pdf(paste0(output_prefix, "_Phylogenetic_Tree_Unrooted.pdf"), width=8, height=6)
    plot(phylo_tree, main="Phylogenetic Reconstruction (Unrooted)", type="unrooted", edge.width=2, cex=1.2)
    dev.off()
})

# --- B. PCA Trajectory ---
cat("LOG: generating Trajectory...\n")
pca <- prcomp(t(mat_sig))
pcaData <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Sample=rownames(meta), Class=meta$Classification)
centroids <- aggregate(cbind(PC1, PC2) ~ Class, data=pcaData, mean)

p_pca <- ggplot(pcaData, aes(x=PC1, y=PC2, color=Class)) +
  geom_point(size=5) +
  geom_path(data=centroids, aes(group=1), arrow=arrow(length=unit(0.3,"cm"), type="closed"), linetype="dashed", color="black", alpha=0.5) +
  geom_text_repel(aes(label=Sample)) +
  scale_color_manual(values=c("blue", "orange", "red")) + theme_bw()
ggsave(paste0(output_prefix, "_PCA_Trajectory.pdf"), p_pca)

# ==============================================================================
# 3. GLOBAL PATHWAY ANALYSIS (GSVA)
# ==============================================================================
cat("LOG: Running Global GSVA...\n")
gmt <- read.gmt(gmt_file)
gene_sets <- split(gmt$gene, gmt$term)
gsva_res <- gsva(mat_vst, gene_sets, method="gsva", kcdf="Gaussian", verbose=FALSE)

path_var <- apply(gsva_res, 1, var)
top_path <- names(sort(path_var, decreasing=TRUE))[1:30]
pdf(paste0(output_prefix, "_GSVA_Global_Heatmap.pdf"), width=10, height=8)
Heatmap(gsva_res[top_path, ], name="Score", column_title="Global Pathway Activity")
dev.off()

# ==============================================================================
# 4. PAIRWISE KITCHEN SINK (Loop through Contrasts)
# ==============================================================================
cat("LOG: Running Pairwise Analysis (Volcano + GSEA)...\n")
dds_wald <- DESeq(dds, test="Wald")
contrasts <- read.csv(contrast_file)

for(i in 1:nrow(contrasts)) {
    ref <- contrasts$reference[i]
    target <- contrasts$target[i]
    comp <- contrasts$id[i]

    cat(sprintf("   -> Contrast: %s vs %s\n", target, ref))
    res <- results(dds_wald, contrast=c("Classification", target, ref))
    res_df <- as.data.frame(res)
    res_df$symbol <- rownames(res_df)

    # --- Enhanced Volcano ---
    p_vol <- EnhancedVolcano(
        res_df,
        lab = res_df$symbol,
        x = 'log2FoldChange',
        y = 'padj',
        title = paste0("Volcano: ", comp),
        pCutoff = 0.05,
        FCcutoff = 1.0,
        pointSize = 2.0,
        labSize = 3.0,
        drawConnectors = TRUE,
        widthConnectors = 0.5,
        max.overlaps = 20
    )
    ggsave(paste0(output_prefix, "_Volcano_", comp, ".pdf"), p_vol, width=8, height=8)

    # --- GSEA ---
    res_df$rank <- sign(res_df$log2FoldChange) * -log10(res_df$pvalue)
    res_df <- res_df[!is.na(res_df$rank) & !is.infinite(res_df$rank),]
    gene_list <- sort(setNames(res_df$rank, res_df$symbol), decreasing=TRUE)

    # Added eps=1e-50 for stability
    gsea_out <- GSEA(gene_list, TERM2GENE=gmt, pvalueCutoff=1, verbose=FALSE, eps=1e-50)

    if(!is.null(gsea_out) && nrow(gsea_out) > 0) {
        gsea_out <- pairwise_termsim(gsea_out)
        p_dot <- dotplot(gsea_out, showCategory=10, split=".sign") + facet_grid(.~.sign)
        ggsave(paste0(output_prefix, "_GSEA_Dot_", comp, ".pdf"), p_dot, width=10, height=8)
    }
}
cat("LOG: Mega-Script Complete.\n")
