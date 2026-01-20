#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# PCA.R: Publication PCA (Hardcoded Positions, Clean Labels, Correct Legends)
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(dplyr)
    library(ggrepel)
    library(EnsDb.Hsapiens.v86)
    library(GSVA)
    library(data.table)
    library(GSEABase)
    library(tidyr)
    library(clusterProfiler)
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
    if (!grepl("_mqc$", filename_base)) filename_base <- paste0(filename_base, "_mqc")
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

# Run GSVA
gsva_sub <- gsva(mat_vst, gbm_sigs, method="gsva", kcdf="Gaussian", verbose=FALSE)

# Normalize scores to percentages for labelling
gsva_norm <- apply(gsva_sub, 2, function(x) {
    x_pos <- x - min(x) + 0.1
    return(round((x_pos / sum(x_pos)) * 100, 1))
})

# ==============================================================================
# 2. PCA WITH HARDCODED LABEL POSITIONS
# ==============================================================================
pca <- prcomp(t(mat_sig))

pcaData <- data.frame(
    PC1=pca$x[,1], 
    PC2=pca$x[,2], 
    PC3=pca$x[,3], 
    Sample=rownames(meta), 
    Class=meta$Classification
)
pcaData$Dominant_Subtype <- apply(gsva_sub, 2, function(x) rownames(gsva_sub)[which.max(x)])

# A. Calculate Class Centroids (Keep for arrows)
centroids <- aggregate(cbind(PC1, PC2) ~ Class, data=pcaData, FUN=mean)
centroids <- centroids[order(centroids$Class),]

# B. Hardcoded Label Coordinates (User Requested)
centroids$Label_X <- NA
centroids$Label_Y <- NA

# Culture: -20, 0
centroids[centroids$Class == "Culture_U2", "Label_X"] <- -20
centroids[centroids$Class == "Culture_U2", "Label_Y"] <- 0

# Primary: -4, 0
centroids[centroids$Class == "Primary_U2", "Label_X"] <- -4
centroids[centroids$Class == "Primary_U2", "Label_Y"] <- 0

# Recurrent: -4, 10
centroids[centroids$Class == "Recurrent_U2", "Label_X"] <- -4
centroids[centroids$Class == "Recurrent_U2", "Label_Y"] <- 10

# C. Generate Rich Labels (Clean Name + Subtype %)
avg_comp <- as.data.frame(t(gsva_norm))
avg_comp$Class <- meta[rownames(avg_comp), "Classification"]
comp_summ <- avg_comp %>% group_by(Class) %>% summarise(across(everything(), mean))

centroids$Label <- apply(centroids, 1, function(row) {
    cls <- row["Class"]
    # CLEAN NAME: Remove _U2 suffix
    clean_name <- sub("_U2", "", cls)
    
    vals <- comp_summ[comp_summ$Class == cls, -1]
    pcts <- paste(names(vals), paste0(round(as.numeric(vals),0), "%"), sep=":", collapse="\n")
    paste0(clean_name, "\n", pcts)
})

arrow_data <- data.frame(x_start = centroids$PC1[-nrow(centroids)], y_start = centroids$PC2[-nrow(centroids)],
                         x_end = centroids$PC1[-1], y_end = centroids$PC2[-1])

# D. Publication Colors
subtype_colors <- c("Mesenchymal"="#E41A1C", "Proneural"="#377EB8", "Classical"="#4DAF4A", "Neural"="#984EA3", "Other/Unclassified"="grey70")

p_pca <- ggplot(pcaData, aes(x=PC1, y=PC2)) +
    # Trajectory Arrow (Connecting actual centroids, not labels)
    geom_segment(data=arrow_data, aes(x=x_start, y=y_start, xend=x_end, yend=y_end),
                 arrow=arrow(length=unit(0.5,"cm"), type="closed"), color="grey60", linewidth=2, inherit.aes=FALSE) +
    
    # Points
    geom_point(aes(fill=Dominant_Subtype, shape=Class), size=9, alpha=0.9, color="black", stroke=1) +
    
    # Labels at Hardcoded Positions
    geom_label(data=centroids, aes(x=Label_X, y=Label_Y, label=Label),
               color="black", fill="white", alpha=0.9, size=4, fontface="bold", inherit.aes=FALSE) +
    
    scale_fill_manual(values=subtype_colors, name="Subtype") +
    scale_shape_manual(values=c(21, 24, 22), labels=function(x) sub("_U2", "", x)) +
    
    # FIX: LEGEND OVERRIDES
    guides(
        # For 'shape' (Class), override fill to be black so shapes are visible
        shape = guide_legend(override.aes = list(fill = "black")),
        # For 'fill' (Subtype), override shape to be 21 (filled circle) so colors show
        fill = guide_legend(override.aes = list(shape = 21, size = 6))
    ) +

    # Clean Theme
    theme_bw(base_size = 18) + 
    theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(vjust = 3),
        axis.title.x = element_text(vjust = -1),
        legend.position = "right",
        # INCREASED VERTICAL SPACING FOR LEGEND KEYS
        legend.key.height = unit(1.2, "cm"),
        plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(title="Evolutionary Trajectory & Subtype Composition", x="PC1", y="PC2")

save_plot(p_pca, paste0(output_prefix, "_PCA_Publication"), 10, 8)

# ==============================================================================
# 3. GEMINI PROMPT GENERATION
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
