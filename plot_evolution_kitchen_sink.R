#!/usr/bin/env Rscript

# ==============================================================================
# U251 EVOLUTION & KITCHEN SINK VISUALIZATION SUITE (MEGA EDITION v3.0)
# ------------------------------------------------------------------------------
# Features:
# 1. ROBUST ID MAPPING: Uses EnsDb.Hsapiens.v86 (Standardized)
# 2. EVOLUTION PLOTS: Phylogenetic Trees (Rooted/Unrooted) & PCA Trajectories
# 3. KITCHEN SINK PLOTS: EnhancedVolcano, GSEA (Dot/Net), GSVA Heatmap
# 4. MULTIQC READY: Saves all outputs as both PDF (vector) and PNG (report)
# ==============================================================================

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(dplyr)
    library(RColorBrewer)
    library(ape)            
    library(ggrepel)      
    library(EnsDb.Hsapiens.v86) # Standardized ID Mapping
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
# 1. DATA LOADING & ROBUST ID MAPPING
# ==============================================================================
cat("LOG: Loading Data...\n")
meta <- read.csv(meta_file, row.names=1)
# Enforce strict biological order
meta$Classification <- factor(meta$Classification, levels = c("Culture_U2", "Primary_U2", "Recurrent_U2"))

counts <- read.table(counts_file, header=TRUE, row.names=1, check.names=FALSE)
counts <- counts[, rownames(meta)]

cat("LOG: Mapping IDs to Symbols (EnsDb)...\n")
clean_ids <- sub("\\..*", "", rownames(counts))
syms <- suppressWarnings(mapIds(EnsDb.Hsapiens.v86,
                                keys=clean_ids,
                                column="SYMBOL",
                                keytype="GENEID",
                                multiVals="first"))
syms[is.na(syms)] <- clean_ids[is.na(syms)]
rownames(counts) <- make.unique(as.character(syms))

dds <- DESeqDataSetFromMatrix(round(counts), meta, ~ Classification)

# ==============================================================================
# 2. EVOLUTIONARY ANALYSIS (LRT) - TREE & TRAJECTORY
# ==============================================================================
cat("LOG: Running Evolution LRT...\n")
dds_lrt <- DESeq(dds, test="LRT", reduced=~1)
res_lrt <- results(dds_lrt)
vsd <- vst(dds_lrt, blind=FALSE)
mat_vst <- assay(vsd)

# --- A. Phylogenetic Tree ---
cat("LOG: Generating Phylogenetic Tree...\n")
sig_genes <- rownames(res_lrt[which(res_lrt$padj < 0.05), ])
mat_sig <- mat_vst[sig_genes, ]
dists <- dist(t(mat_sig))
phylo_tree <- as.phylo(hclust(dists, method="ward.D2"))

tryCatch({
    # Attempt Rooting at Culture (C2B)
    rooted_tree <- root(phylo_tree, outgroup="C2B", resolve.root=TRUE)
    
    # Save PDF & PNG
    pdf(paste0(output_prefix, "_Phylogenetic_Tree_mqc.pdf"), width=8, height=6)
    plot(rooted_tree, main="Phylogenetic Reconstruction", type="phylogram", edge.width=2, cex=1.2)
    tiplabels(pch=21, bg=c("blue", rep("orange", 3), rep("red", 3)), cex=2)
    dev.off()
    
    png(paste0(output_prefix, "_Phylogenetic_Tree_mqc.png"), width=8, height=6, units="in", res=300)
    plot(rooted_tree, main="Phylogenetic Reconstruction", type="phylogram", edge.width=2, cex=1.2)
    tiplabels(pch=21, bg=c("blue", rep("orange", 3), rep("red", 3)), cex=2)
    dev.off()
}, error = function(e) {
    cat("WARN: Could not root tree at C2B. Plotting unrooted.\n")
    # Fallback to Unrooted
    pdf(paste0(output_prefix, "_Phylogenetic_Tree_Unrooted_mqc.pdf"), width=8, height=6)
    plot(phylo_tree, type="unrooted")
    dev.off()
    
    png(paste0(output_prefix, "_Phylogenetic_Tree_Unrooted_mqc.png"), width=8, height=6, units="in", res=300)
    plot(phylo_tree, type="unrooted")
    dev.off()
})

# --- B. PCA Trajectory ---
cat("LOG: Generating PCA Trajectory...\n")
pca <- prcomp(t(mat_sig))
pcaData <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Sample=rownames(meta), Class=meta$Classification)
centroids <- aggregate(cbind(PC1, PC2) ~ Class, data=pcaData, mean)

p_pca <- ggplot(pcaData, aes(x=PC1, y=PC2, color=Class)) +
  geom_point(size=5) +
  geom_path(data=centroids, aes(group=1), arrow=arrow(length=unit(0.3,"cm"), type="closed"), linetype="dashed", color="black", alpha=0.5) +
  geom_text_repel(aes(label=Sample)) +
  scale_color_manual(values=c("blue", "orange", "red")) + theme_bw()

ggsave(paste0(output_prefix, "_PCA_Trajectory_mqc.pdf"), p_pca, width=8, height=6)
ggsave(paste0(output_prefix, "_PCA_Trajectory_mqc.png"), p_pca, width=8, height=6, dpi=300)

# ==============================================================================
# 3. GLOBAL PATHWAY ANALYSIS (GSVA)
# ==============================================================================
cat("LOG: Running Global GSVA...\n")
gmt <- read.gmt(gmt_file)
gene_sets <- split(gmt$gene, gmt$term)
gsva_res <- gsva(mat_vst, gene_sets, method="gsva", kcdf="Gaussian", verbose=FALSE)

path_var <- apply(gsva_res, 1, var)
top_path <- names(sort(path_var, decreasing=TRUE))[1:30]
ht_global <- Heatmap(gsva_res[top_path, ], name="Score", column_title="Global Pathway Activity")

pdf(paste0(output_prefix, "_GSVA_Global_Heatmap_mqc.pdf"), width=10, height=8)
draw(ht_global)
dev.off()

png(paste0(output_prefix, "_GSVA_Global_Heatmap_mqc.png"), width=10, height=8, units="in", res=300)
draw(ht_global)
dev.off()

# ==============================================================================
# 4. PAIRWISE CONTRASTS (VOLCANO + GSEA)
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
    ggsave(paste0(output_prefix, "_Volcano_", comp, "_mqc.pdf"), p_vol, width=8, height=8)
    ggsave(paste0(output_prefix, "_Volcano_", comp, "_mqc.png"), p_vol, width=8, height=8, dpi=300)
    
    # --- GSEA ---
    res_df$rank <- sign(res_df$log2FoldChange) * -log10(res_df$pvalue)
    res_df <- res_df[!is.na(res_df$rank) & !is.infinite(res_df$rank),]
    gene_list <- sort(setNames(res_df$rank, res_df$symbol), decreasing=TRUE)
    
    # Run GSEA with stability fix
    gsea_out <- GSEA(gene_list, TERM2GENE=gmt, pvalueCutoff=1, verbose=FALSE, eps=1e-50)
    
    if(!is.null(gsea_out) && nrow(gsea_out) > 0) {
        if(nrow(gsea_out) > 50) gsea_out@result <- head(gsea_out@result[order(gsea_out@result$p.adjust), ], 50)
        gsea_out <- pairwise_termsim(gsea_out)
        
        # Dotplot
        p_dot <- dotplot(gsea_out, showCategory=10, split=".sign") + facet_grid(.~.sign)
        ggsave(paste0(output_prefix, "_GSEA_Dot_", comp, "_mqc.pdf"), p_dot, width=10, height=8)
        ggsave(paste0(output_prefix, "_GSEA_Dot_", comp, "_mqc.png"), p_dot, width=10, height=8, dpi=300)
        
        # Cnet Plot (Network)
        p_cnet <- cnetplot(gsea_out, categorySize="pvalue", foldChange=gene_list, circular=TRUE, color.params=list(edge=TRUE))
        ggsave(paste0(output_prefix, "_GSEA_Network_", comp, "_mqc.pdf"), p_cnet, width=12, height=12)
        ggsave(paste0(output_prefix, "_GSEA_Network_", comp, "_mqc.png"), p_cnet, width=12, height=12, dpi=300)
    }
}
cat("LOG: Mega-Script Complete.\n")
