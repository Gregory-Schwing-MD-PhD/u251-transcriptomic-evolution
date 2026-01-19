#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# PCA.R: Publication PCA (Nudged Labels), 3D Trajectory, & AI Prompts
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(dplyr)
    library(ggrepel)
    library(EnsDb.Hsapiens.v86)
    library(GSVA)
    library(scatterplot3d) 
    library(data.table)
    library(GSEABase)
    library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: PCA.R <counts> <meta> <gmt> <out_prefix>")

counts_file   <- args[1]
meta_file     <- args[2]
gmt_file      <- args[3]
output_prefix <- args[4]

# --- Setup ---
dir.create(dirname(output_prefix), showWarnings = FALSE, recursive = TRUE)

save_plot <- function(plot_obj, filename_base, w=10, h=8) {
    pdf(paste0(filename_base, ".pdf"), width=w, height=h)
    print(plot_obj)
    dev.off()
    png(paste0(filename_base, ".png"), width=w, height=h, units="in", res=300)
    print(plot_obj)
    dev.off()
}

# --- Load Data ---
meta <- read.csv(meta_file, row.names=1)
if("Classification" %in% colnames(meta)) {
    meta$Classification <- factor(meta$Classification, levels = c("Culture_U2", "Primary_U2", "Recurrent_U2"))
}

counts <- read.table(counts_file, header=TRUE, row.names=1, check.names=FALSE)
counts <- counts[, rownames(meta)]
clean_ids <- sub("\\..*", "", rownames(counts))
syms <- suppressWarnings(mapIds(EnsDb.Hsapiens.v86, keys=clean_ids, column="SYMBOL", keytype="GENEID", multiVals="first"))
syms[is.na(syms)] <- clean_ids[is.na(syms)]
rownames(counts) <- make.unique(as.character(syms))

# --- VST Transform ---
dds <- DESeqDataSetFromMatrix(round(counts), meta, ~1)
vsd <- vst(dds, blind=FALSE)
mat_vst <- assay(vsd)
mat_sig <- mat_vst[head(order(rowVars(mat_vst), decreasing=TRUE), 500), ]

# ==============================================================================
# 1. GBM SUBTYPING (GSVA)
# ==============================================================================
gbm_sigs <- list(
    Mesenchymal = c("CHI3L1", "CD44", "VIM", "RELB", "STAT3", "MET", "TRADD", "MMP9", "TIMP1"),
    Classical   = c("EGFR", "AKT2", "NOTCH3", "JAG1", "CCND2", "F3", "PDGFA", "NES"),
    Proneural   = c("PDGFRA", "IDH1", "OLIG2", "SOX2", "NKX2-2", "OLIG1", "TP53"),
    Neural      = c("NEFL", "GABRA1", "SYT1", "SLC12A5", "MBP", "GABRG2")
)
gsva_sub <- gsva(mat_vst, gbm_sigs, method="gsva", kcdf="Gaussian", verbose=FALSE)

# ==============================================================================
# 2. PCA WITH CENTER-NUDGED LABELS (PUBLICATION STYLE)
# ==============================================================================
pca <- prcomp(t(mat_sig))
pcaData <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Sample=rownames(meta), Class=meta$Classification)
pcaData$Dominant_Subtype <- apply(gsva_sub, 2, function(x) rownames(gsva_sub)[which.max(x)])

# Calculate Global Center & Centroids
global_center_x <- mean(pcaData$PC1)
global_center_y <- mean(pcaData$PC2)
centroids <- aggregate(cbind(PC1, PC2) ~ Class, data=pcaData, FUN=mean)

# Nudge labels 50% towards global center
nudge_factor <- 0.5 
centroids$Label_X <- centroids$PC1 + (global_center_x - centroids$PC1) * nudge_factor
centroids$Label_Y <- centroids$PC2 + (global_center_y - centroids$PC2) * nudge_factor

arrow_data <- data.frame(x_start = centroids$PC1[-nrow(centroids)], y_start = centroids$PC2[-nrow(centroids)],
                         x_end = centroids$PC1[-1], y_end = centroids$PC2[-1])

# Colors for Publication (Distinct High Contrast)
pub_colors <- c("Culture_U2"="#1B9E77", "Primary_U2"="#D95F02", "Recurrent_U2"="#7570B3")
subtype_colors <- c("Mesenchymal"="#E41A1C", "Proneural"="#377EB8", "Classical"="#4DAF4A", "Neural"="#984EA3")

p_pca <- ggplot(pcaData, aes(x=PC1, y=PC2)) +
    geom_segment(data=arrow_data, aes(x=x_start, y=y_start, xend=x_end, yend=y_end),
                 arrow=arrow(length=unit(0.5,"cm"), type="closed"), color="grey60", size=2, inherit.aes=FALSE) +
    
    # Points with strokes for clarity
    geom_point(aes(fill=Dominant_Subtype, shape=Class), size=9, alpha=0.9, color="black", stroke=1) +
    
    # Nudged Centroid Labels
    geom_label(data=centroids, aes(x=Label_X, y=Label_Y, label=Class),
               color="black", fill="white", alpha=0.9, size=6, fontface="bold", inherit.aes=FALSE) +
    
    # Manual Scales
    scale_fill_manual(values=subtype_colors, name="Subtype") +
    scale_shape_manual(values=c(21, 24, 22)) +
    
    # Publication Theme (Clean, No Grid)
    theme_bw(base_size = 20) + 
    theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(vjust = 3),
        axis.title.x = element_text(vjust = -1),
        legend.position = "right",
        plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(title="Evolutionary Trajectory", x="PC1", y="PC2")

save_plot(p_pca, paste0(output_prefix, "_PCA_Publication"), 10, 8)

# ==============================================================================
# 3. 3D TRAJECTORY
# ==============================================================================
data_3d <- as.data.frame(t(gsva_sub[c("Classical", "Neural", "Mesenchymal"), ]))
data_3d$Class <- meta[rownames(data_3d), "Classification"]
point_colors <- pub_colors[as.character(data_3d$Class)]

pdf(paste0(output_prefix, "_3D_Subtype_Trajectory.pdf"), width=10, height=10)
s3d <- scatterplot3d(data_3d$Classical, data_3d$Neural, data_3d$Mesenchymal,
              pch=16, color=point_colors, size=2, type="h",
              main="3D Subtype Evolution",
              xlab="Classical Score", ylab="Neural Score", zlab="Mesenchymal Score",
              grid=TRUE, box=FALSE) # Removed box for cleaner look
legend(s3d$xyz.convert(max(data_3d$Classical), max(data_3d$Neural), min(data_3d$Mesenchymal)), 
       legend = names(pub_colors), col = pub_colors, pch = 16, cex=1.5, bty="n")
centroids_3d <- aggregate(cbind(Classical, Neural, Mesenchymal) ~ Class, data=data_3d, FUN=mean)
s3d$points3d(centroids_3d$Classical, centroids_3d$Neural, centroids_3d$Mesenchymal, 
             col="black", type="l", lwd=3, lty=2)
dev.off()

png(paste0(output_prefix, "_3D_Subtype_Trajectory.png"), width=10, height=10, units="in", res=300)
s3d <- scatterplot3d(data_3d$Classical, data_3d$Neural, data_3d$Mesenchymal,
              pch=16, color=point_colors, size=2, type="h",
              main="3D Subtype Evolution",
              xlab="Classical Score", ylab="Neural Score", zlab="Mesenchymal Score",
              grid=TRUE, box=FALSE)
legend(s3d$xyz.convert(max(data_3d$Classical), max(data_3d$Neural), min(data_3d$Mesenchymal)), 
       legend = names(pub_colors), col = pub_colors, pch = 16, cex=1.5, bty="n")
s3d$points3d(centroids_3d$Classical, centroids_3d$Neural, centroids_3d$Mesenchymal, 
             col="black", type="l", lwd=3, lty=2)
dev.off()

# ==============================================================================
# 4. GEMINI PROMPT GENERATION
# ==============================================================================
loadings <- as.data.frame(pca$rotation[, 1:5])
top_drivers <- list()
for(i in 1:5) {
    pc <- paste0("PC", i)
    genes <- rownames(loadings[order(abs(loadings[[pc]]), decreasing=TRUE), ][1:10, ])
    top_drivers[[pc]] <- paste(genes, collapse=", ")
}

gmt <- read.gmt(gmt_file)
gsva_global <- gsva(mat_vst, split(gmt$gene, gmt$term), method="gsva", kcdf="Gaussian", verbose=FALSE)

sample_summaries <- list()
for(s in rownames(meta)) {
    cls <- meta[s, "Classification"]
    top_paths <- names(sort(gsva_global[,s], decreasing=TRUE)[1:3])
    pc_coords <- paste(round(pcaData[s, c("PC1", "PC2", "PC3")], 2), collapse=", ")
    sample_summaries[[s]] <- paste0("Sample: ", s, " (", cls, ") | Top Paths: ", paste(top_paths, collapse=", "), " | PCA: ", pc_coords)
}

prompt_file <- paste0(dirname(output_prefix), "/GEMINI_PROMPTS.txt")
sink(prompt_file)
cat("I have transcriptomic data for U251 GBM evolution (Culture -> Primary -> Recurrent).\n")
cat("Please interpret the biological trajectory based on these PCA drivers and Sample Phenotypes.\n\n")
cat("=== PART 1: TOP PCA GENES ===\n")
for(pc in names(top_drivers)) cat(paste0("- ", pc, ": ", top_drivers[[pc]], "\n"))
cat("\n=== PART 2: SAMPLE PROFILES ===\n")
cat(paste(unlist(sample_summaries), collapse="\n"))
sink()
