#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# EVOLUTIONARY KITCHEN SINK v7.0 (With GBM Subtyping)
# Combines: Evolution (Tree/PCA), GBM Subtyping, GSEA, GSVA, Volcano
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(dplyr)
    library(RColorBrewer)
    library(ape)
    library(ggrepel)
    library(EnsDb.Hsapiens.v86)
    library(clusterProfiler)
    library(enrichplot)
    library(GSVA)
    library(GSEABase)
    library(ComplexHeatmap)
    library(circlize)
    library(EnhancedVolcano)
    library(ggnewscale)
    library(data.table)
    if (requireNamespace("ggupset", quietly = TRUE)) library(ggupset)
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

# Helper for saving
save_plot <- function(plot_obj, filename_base, w=10, h=8) {
    tryCatch({
        if (inherits(plot_obj, "Heatmap")) {
            pdf(paste0(filename_base, ".pdf"), width=w, height=h); draw(plot_obj); dev.off()
            png(paste0(filename_base, ".png"), width=w, height=h, units="in", res=300); draw(plot_obj); dev.off()
        } else {
            ggsave(paste0(filename_base, ".pdf"), plot_obj, width=w, height=h)
            ggsave(paste0(filename_base, ".png"), plot_obj, width=w, height=h, dpi=300, bg="white")
        }
    }, error = function(e) cat(paste("WARN: Save failed for", filename_base, "\n")))
}

# ==============================================================================
# 1. DATA LOADING & MAPPING
# ==============================================================================
cat("LOG [1/7]: Loading Data & Mapping IDs...\n")
meta <- read.csv(meta_file, row.names=1)

# Ensure factors (Culture -> Primary -> Recurrent)
if("Classification" %in% colnames(meta)) {
    meta$Classification <- factor(meta$Classification, levels = c("Culture_U2", "Primary_U2", "Recurrent_U2"))
}

counts <- read.table(counts_file, header=TRUE, row.names=1, check.names=FALSE)
counts <- counts[, rownames(meta)]

# Map IDs to Symbols
clean_ids <- sub("\\..*", "", rownames(counts))
syms <- suppressWarnings(mapIds(EnsDb.Hsapiens.v86, keys=clean_ids, column="SYMBOL", keytype="GENEID", multiVals="first"))
syms[is.na(syms)] <- clean_ids[is.na(syms)]
rownames(counts) <- make.unique(as.character(syms))

# Design
contrasts <- read.csv(contrast_file)
design_var <- unique(contrasts$variable)[1]
design_formula <- as.formula(paste0("~ ", design_var))
dds <- DESeqDataSetFromMatrix(round(counts), meta, design_formula)

# ==============================================================================
# 2. EVOLUTIONARY RECONSTRUCTION
# ==============================================================================
cat("LOG [2/7]: Evolutionary Reconstruction (Tree + PCA)...\n")
dds_lrt <- DESeq(dds, test="LRT", reduced=~1)
res_lrt <- results(dds_lrt)
vsd <- vst(dds_lrt, blind=FALSE)
mat_vst <- assay(vsd)

# Colors
cls_colors <- c("Culture_U2"="blue", "Primary_U2"="orange", "Recurrent_U2"="red")
if (!all(levels(meta$Classification) %in% names(cls_colors))) {
  cls_colors <- setNames(rainbow(length(levels(meta$Classification))), levels(meta$Classification))
}

# --- A. Phylogenetic Tree ---
sig_genes <- rownames(res_lrt[which(res_lrt$padj < 0.05), ])
if(length(sig_genes) < 5) sig_genes <- head(order(rowVars(mat_vst), decreasing=TRUE), 500)
mat_sig <- mat_vst[sig_genes, ]
phylo_tree <- as.phylo(hclust(dist(t(mat_sig)), method="ward.D2"))

tryCatch({
    if("C2B" %in% phylo_tree$tip.label) rooted_tree <- root(phylo_tree, outgroup="C2B", resolve.root=TRUE) else rooted_tree <- phylo_tree
    tip_cols <- cls_colors[as.character(meta[rooted_tree$tip.label, "Classification"])]
    
    pdf(paste0(output_prefix, "_Phylogenetic_Tree_mqc.pdf"), width=8, height=6)
    plot(rooted_tree, main="Phylogenetic Reconstruction", type="phylogram", edge.width=2, tip.color=tip_cols, label.offset=0.5)
    legend("topleft", legend=names(cls_colors), fill=cls_colors, bty="n", cex=0.8)
    dev.off()
    
    png(paste0(output_prefix, "_Phylogenetic_Tree_mqc.png"), width=8, height=6, units="in", res=300)
    plot(rooted_tree, main="Phylogenetic Reconstruction", type="phylogram", edge.width=2, tip.color=tip_cols, label.offset=0.5)
    legend("topleft", legend=names(cls_colors), fill=cls_colors, bty="n", cex=0.8)
    dev.off()
}, error = function(e) cat("WARN: Tree failed.\n"))

# --- B. PCA Trajectory ---
pca <- prcomp(t(mat_sig))
pcaData <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Sample=rownames(meta), Class=meta$Classification)
centroids <- aggregate(cbind(PC1, PC2) ~ Class, data=pcaData, FUN=mean)
centroids <- centroids[order(centroids$Class),] 

arrow_data <- data.frame(x_start = centroids$PC1[-nrow(centroids)], y_start = centroids$PC2[-nrow(centroids)],
                         x_end = centroids$PC1[-1], y_end = centroids$PC2[-1])

p_pca <- ggplot(pcaData, aes(x=PC1, y=PC2, color=Class)) +
  geom_point(size=5, alpha=0.8) +
  geom_text_repel(aes(label=Sample), show.legend=FALSE) +
  geom_segment(data=arrow_data, aes(x=x_start, y=y_start, xend=x_end, yend=y_end), 
               arrow=arrow(length=unit(0.3,"cm"), type="closed"), color="black", size=1, inherit.aes=FALSE) +
  scale_color_manual(values=cls_colors) + theme_bw() + labs(title="PCA Trajectory")
save_plot(p_pca, paste0(output_prefix, "_PCA_Trajectory_mqc"), 8, 6)

# ==============================================================================
# 3. GBM SUBTYPE CLASSIFICATION (NEW MODULE) 
# ==============================================================================
cat("LOG [3/7]: Running GBM Subtype Classification (Verhaak/Wang)...\n")

# Hardcoded Canonical Signatures (Robust fallback)
gbm_sigs <- list(
    Mesenchymal = c("CHI3L1", "CD44", "VIM", "RELB", "STAT3", "MET", "TRADD", "MMP9", "TIMP1"),
    Classical   = c("EGFR", "AKT2", "NOTCH3", "JAG1", "CCND2", "F3", "PDGFA", "NES"),
    Proneural   = c("PDGFRA", "IDH1", "OLIG2", "SOX2", "NKX2-2", "OLIG1", "TP53"),
    Neural      = c("NEFL", "GABRA1", "SYT1", "SLC12A5", "MBP", "GABRG2")
)

# Run GSVA for Subtypes
gsva_sub <- gsva(mat_vst, gbm_sigs, method="gsva", kcdf="Gaussian", verbose=FALSE)

# --- Plot 1: Subtype Heatmap ---
ht_sub <- Heatmap(gsva_sub, name="Subtype Score", column_title="GBM Subtype Classification",
                  col=colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")),
                  top_annotation = HeatmapAnnotation(Class = meta$Classification, 
                                                     col = list(Class = cls_colors)))
save_plot(ht_sub, paste0(output_prefix, "_Subtype_Heatmap_mqc"), 8, 5)

# --- Plot 2: Subtype Evolution Line Plot ---
# Reshape for ggplot
sub_df <- as.data.frame(t(gsva_sub))
sub_df$Sample <- rownames(sub_df)
sub_df$Class <- meta[sub_df$Sample, "Classification"]
sub_long <- reshape2::melt(sub_df, id.vars=c("Sample", "Class"), variable.name="Subtype", value.name="Score")

# Calculate Mean Score per Class per Subtype
sub_summ <- sub_long %>% group_by(Class, Subtype) %>% summarize(MeanScore = mean(Score), .groups="drop")

p_evol <- ggplot(sub_summ, aes(x=Class, y=MeanScore, color=Subtype, group=Subtype)) +
    geom_line(size=1.5) + geom_point(size=4) +
    scale_color_brewer(palette="Set1") +
    theme_bw() +
    labs(title="Evolution of GBM Subtypes", y="Mean GSVA Score", x="Evolutionary Stage") +
    theme(axis.text.x = element_text(angle=45, hjust=1))

save_plot(p_evol, paste0(output_prefix, "_Subtype_Evolution_mqc"), 8, 6)

# ==============================================================================
# 4. GLOBAL PATHWAYS (GSVA)
# ==============================================================================
cat("LOG [4/7]: Global GSVA...\n")
gmt <- read.gmt(gmt_file)
gsva_res <- gsva(mat_vst, split(gmt$gene, gmt$term), method="gsva", kcdf="Gaussian", verbose=FALSE)
top_path <- names(sort(apply(gsva_res, 1, var), decreasing=TRUE)[1:30])
ht_global <- Heatmap(gsva_res[top_path, ], name="Score", column_title="Global Pathway Activity",
                     col=colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")))
save_plot(ht_global, paste0(output_prefix, "_GSVA_Global_Heatmap_mqc"), 10, 8)

# ==============================================================================
# 5. PAIRWISE CONTRASTS
# ==============================================================================
cat("LOG [5/7]: Pairwise Analysis...\n")
dds_wald <- DESeq(dds, test="Wald")

for(i in 1:nrow(contrasts)) {
    comp <- contrasts$id[i]; var <- contrasts$variable[i]; ref <- contrasts$reference[i]; target <- contrasts$target[i]
    cat(sprintf("   -> Contrast: %s\n", comp))
    
    res <- results(dds_wald, contrast=c(var, target, ref))
    res_df <- as.data.frame(res); res_df$symbol <- rownames(res_df)
    
    # Volcano
    p_vol <- EnhancedVolcano(res_df, lab=res_df$symbol, x='log2FoldChange', y='pvalue',
        title=paste0("Volcano: ", comp), subtitle=paste0(target, " vs ", ref),
        pCutoff=0.05, FCcutoff=1.0, pointSize=2.0, labSize=3.0, drawConnectors=TRUE, widthConnectors=0.5, max.overlaps=20)
    save_plot(p_vol, paste0(output_prefix, "_Volcano_", comp, "_mqc"), 8, 8)
    
    # GSEA
    res_df$rank <- sign(res_df$log2FoldChange) * -log10(res_df$pvalue)
    res_df <- res_df[!is.na(res_df$rank) & !is.infinite(res_df$rank),]
    gene_list <- sort(setNames(res_df$rank, res_df$symbol), decreasing=TRUE)
    
    gsea_out <- GSEA(gene_list, TERM2GENE=gmt, pvalueCutoff=1, verbose=FALSE, eps=1e-50)
    if(!is.null(gsea_out) && nrow(gsea_out) > 0) {
        if(nrow(gsea_out) > 50) gsea_out@result <- head(gsea_out@result[order(gsea_out@result$p.adjust), ], 50)
        gsea_out <- pairwise_termsim(gsea_out)
        
        save_plot(dotplot(gsea_out, showCategory=15, split=".sign") + facet_grid(.~.sign), paste0(output_prefix, "_GSEA_Dot_", comp, "_mqc"), 10, 8)
        tryCatch({ save_plot(cnetplot(gsea_out, categorySize="pvalue", foldChange=gene_list, circular=TRUE), paste0(output_prefix, "_GSEA_Cnet_", comp, "_mqc"), 12, 12) }, error=function(e) NULL)
        tryCatch({ save_plot(emapplot(gsea_out, showCategory=20, cex.params=list(category_label=0.6)), paste0(output_prefix, "_GSEA_Emap_", comp, "_mqc"), 12, 10) }, error=function(e) NULL)
        if (requireNamespace("ggupset", quietly=TRUE)) tryCatch({ save_plot(upsetplot(gsea_out, n=10), paste0(output_prefix, "_GSEA_UpSet_", comp, "_mqc"), 12, 6) }, error=function(e) NULL)
    }
}
cat("LOG [6/7]: Complete.\n")
