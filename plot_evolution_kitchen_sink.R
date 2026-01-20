#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# EVOLUTIONARY KITCHEN SINK v15.7 (Constant Subtype Order)
# ------------------------------------------------------------------------------
# 1. PCA/Tree: Uses 'all.vst.tsv' (Pipeline VST) -> Calculates Trajectory & Subtype %
# 2. Volcano:  Uses 'deseq2.results.tsv' (Pipeline DE)
# 3. GSEA:     Uses Log2FC (Pipeline) -> Dotplots (No Cnet)
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(ape); library(ggrepel)
    library(EnsDb.Hsapiens.v86); library(clusterProfiler); library(enrichplot)
    library(GSVA); library(GSEABase); library(ComplexHeatmap); library(circlize)
    library(EnhancedVolcano); library(data.table); library(tidyr)
    if (requireNamespace("ggupset", quietly = TRUE)) library(ggupset)
})

# --- ARGS ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Usage: script.R <results_dir> <out_prefix> <meta> <contrasts> <gmt>")

results_dir   <- args[1]
output_prefix <- args[2]
meta_file     <- args[3]
contrast_file <- args[4]
gmt_file      <- args[5]

# Output setup
dir.create(dirname(output_prefix), showWarnings = FALSE, recursive = TRUE)
gemini_dir <- paste0(dirname(output_prefix), "/GEMINI_DATA")
subtype_dir <- paste0(gemini_dir, "/subtypes")
dir.create(gemini_dir, showWarnings = FALSE)
dir.create(subtype_dir, showWarnings = FALSE)

# --- SAVE FUNCTION ---
# CHANGED: Default size to 9x8
save_plot <- function(plot_obj, filename_base, w=9, h=8) {
    tryCatch({
        if (inherits(plot_obj, "Heatmap")) {
            pdf(paste0(filename_base, ".pdf"), width=w, height=h); draw(plot_obj); dev.off()
            png(paste0(filename_base, ".png"), width=w, height=h, units="in", res=300); draw(plot_obj); dev.off()
        } else {
            ggsave(paste0(filename_base, ".pdf"), plot_obj, width=w, height=h)
            ggsave(paste0(filename_base, ".png"), plot_obj, width=w, height=h, dpi=300, bg="white")
        }
        cat(paste0("SUCCESS: Saved ", basename(filename_base), "\n"))
    }, error = function(e) cat(paste0("ERROR: Failed to save ", basename(filename_base), ": ", e$message, "\n")))
}

# ==============================================================================
# 1. LOAD PRE-COMPUTED VST
# ==============================================================================
cat("LOG [1/5]: Loading Pipeline VST Matrix...\n")
vst_file <- file.path(results_dir, "tables/processed_abundance/all.vst.tsv")

if (!file.exists(vst_file)) stop(paste("FATAL: Pipeline VST file not found at:", vst_file))

# Load VST
mat_vst <- as.matrix(read.table(vst_file, header=TRUE, row.names=1, check.names=FALSE))

# Load Metadata & Sync
meta <- read.csv(meta_file, row.names=1)
common <- intersect(colnames(mat_vst), rownames(meta))
mat_vst <- mat_vst[, common]; meta <- meta[common, , drop=FALSE]

# CLEAN LABELS: Remove _U2 suffix
if("Classification" %in% colnames(meta)) {
    # Ensure correct order for trajectory
    meta$Classification <- factor(meta$Classification, levels=c("Culture_U2", "Primary_U2", "Recurrent_U2"))
}

# ==============================================================================
# 2. EVOLUTIONARY TREE (Sample Clustering)
# ==============================================================================
cat("LOG [2/5]: Generating Tree...\n")
top_var <- head(order(apply(mat_vst, 1, var), decreasing=TRUE), 500)
mat_sig <- mat_vst[top_var, ]

phylo_tree <- as.phylo(hclust(dist(t(mat_sig)), method="ward.D2"))
tryCatch({
    if("C2B" %in% phylo_tree$tip.label) rooted_tree <- root(phylo_tree, outgroup="C2B", resolve.root=TRUE) else rooted_tree <- phylo_tree
    # Size updated via save_plot default, but Tree usually calls pdf directly, updating here too:
    pdf(paste0(output_prefix, "_Phylogenetic_Tree_mqc.pdf"), width=9, height=8)
    plot(rooted_tree, main="Evolutionary Reconstruction", type="phylogram", edge.width=2, cex=1.0)
    dev.off()
}, error = function(e) cat("WARN: Tree failed.\n"))

# ==============================================================================
# 3. PCA TRAJECTORY (With Info Boxes)
# ==============================================================================
cat("LOG [3/5]: Generating PCA with Trajectory Boxes...\n")
cat("   -> Mapping Ensembl IDs to Gene Symbols for GSVA...\n")

clean_ids <- sub("\\..*", "", rownames(mat_vst))
mapped_syms <- mapIds(EnsDb.Hsapiens.v86, keys=clean_ids, column="SYMBOL", keytype="GENEID", multiVals="first")

mat_vst_sym <- mat_vst
rownames(mat_vst_sym) <- ifelse(is.na(mapped_syms), rownames(mat_vst), mapped_syms)
mat_vst_sym <- mat_vst_sym[!is.na(rownames(mat_vst_sym)), ]
mat_vst_sym <- mat_vst_sym[!duplicated(rownames(mat_vst_sym)), ]

# Signatures
gbm_sigs <- list(
    Mesenchymal = c("CHI3L1", "CD44", "VIM", "RELB", "STAT3", "MET", "TRADD", "MMP9", "TIMP1"),
    Classical   = c("EGFR", "AKT2", "NOTCH3", "JAG1", "CCND2", "F3", "PDGFA", "NES"),
    Proneural   = c("PDGFRA", "IDH1", "OLIG2", "SOX2", "NKX2-2", "OLIG1", "TP53"),
    Neural      = c("NEFL", "GABRA1", "SYT1", "SLC12A5", "MBP", "GABRG2")
)

# 1. GSVA & Normalization
gsva_sub <- gsva(mat_vst_sym, gbm_sigs, method="gsva", kcdf="Gaussian", verbose=FALSE)

# Normalize to percentages (positive shift)
gsva_norm <- apply(gsva_sub, 2, function(x) {
    x_pos <- x - min(x) + 0.1
    return(round((x_pos / sum(x_pos)) * 100, 0))
})

# 2. PCA Calculation
pca <- prcomp(t(mat_sig))
pcaData <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Sample=rownames(meta), Class=meta$Classification)
pcaData$Dominant_Subtype <- apply(gsva_sub, 2, function(x) rownames(gsva_sub)[which.max(x)])

# 3. Calculate Centroids & Labels
centroids <- aggregate(cbind(PC1, PC2) ~ Class, data=pcaData, FUN=mean)
centroids <- centroids[order(centroids$Class),] # Culture -> Primary -> Recurrent

# Hardcoded Label Positions
centroids$Label_X <- NA; centroids$Label_Y <- NA
centroids[centroids$Class == "Culture_U2", c("Label_X", "Label_Y")] <- c(-20, 5)
centroids[centroids$Class == "Primary_U2", c("Label_X", "Label_Y")] <- c(-5, 5)
centroids[centroids$Class == "Recurrent_U2", c("Label_X", "Label_Y")] <- c(-5, -10)

# Generate Box Text
avg_comp <- as.data.frame(t(gsva_norm))
avg_comp$Class <- meta[rownames(avg_comp), "Classification"]
# Summarize mean % per class
comp_summ <- avg_comp %>% group_by(Class) %>% summarise(across(everything(), mean))

centroids$Label <- apply(centroids, 1, function(row) {
    cls <- row["Class"]
    # CLEAN NAME: Remove _U2 suffix
    clean_name <- sub("_U2", "", cls)
    
    # Get percentages
    vals <- comp_summ[comp_summ$Class == cls, -1] # Exclude Class col
    
    # CHANGED: REMOVED SORTING (Now Constant Order)
    # vals <- vals[, order(as.numeric(vals), decreasing=TRUE)]
    
    pcts <- paste(names(vals), paste0(round(as.numeric(vals),0), "%"), sep=": ", collapse="\n")
    # CHANGED: Removed dashed line, just clean text
    paste0(clean_name, "\n", pcts)
})

# Trajectory Arrows (Centroid to Centroid)
arrow_data <- data.frame(
    x_start = centroids$PC1[-nrow(centroids)], 
    y_start = centroids$PC2[-nrow(centroids)],
    x_end = centroids$PC1[-1], 
    y_end = centroids$PC2[-1]
)

# 4. Plotting
subtype_colors <- c("Mesenchymal"="#E41A1C", "Proneural"="#377EB8", "Classical"="#4DAF4A", "Neural"="#984EA3", "Other/Unclassified"="grey70")

p_pca <- ggplot(pcaData, aes(x=PC1, y=PC2)) +
    # Trajectory Line
    geom_segment(data=arrow_data, aes(x=x_start, y=y_start, xend=x_end, yend=y_end),
                 arrow=arrow(length=unit(0.4,"cm"), type="closed"), color="grey50", linewidth=1.5, inherit.aes=FALSE) +
    
    # Points
    geom_point(aes(fill=Dominant_Subtype, shape=Class), size=7, color="black", stroke=0.8) +
    
    # Info Boxes
    # CHANGED: fontface="plain" to avoid bolding the subtypes
    geom_label(data=centroids, aes(x=Label_X, y=Label_Y, label=Label),
               fill="white", alpha=0.9, size=3.5, fontface="plain", inherit.aes=FALSE) +
    
    scale_fill_manual(values=subtype_colors) +
    scale_shape_manual(values=c(21, 24, 22), labels=function(x) sub("_U2", "", x)) +
    
    # Guides Override for visibility
    guides(
        shape = guide_legend(override.aes = list(fill = "black")),
        fill = guide_legend(override.aes = list(shape = 21, size = 5))
    ) +
    
    theme_bw(base_size=14) +
    theme(
        legend.position="right",
        legend.spacing.y = unit(0.8, 'cm'),
        legend.key.height = unit(1.2, 'cm'), # ADDED: Fixes icon overlap
        panel.grid.minor = element_blank()
    ) +
    labs(title="Evolutionary Trajectory", subtitle="Pipeline VST | Class Composition")

save_plot(p_pca, paste0(output_prefix, "_PCA_Publication_mqc"), 9, 8)

# ==============================================================================
# 4. CONTRAST VISUALIZATION
# ==============================================================================
cat("LOG [4/5]: Visualizing Pipeline DE Tables...\n")
contrasts <- read.csv(contrast_file)
gmt <- read.gmt(gmt_file)
prompt_summary <- list()

for(i in 1:nrow(contrasts)) {
    cid <- contrasts$id[i]; target <- contrasts$target[i]; ref <- contrasts$reference[i]
    res_file <- file.path(results_dir, "tables/differential", paste0(cid, ".deseq2.results.tsv"))
    
    cat(sprintf("    -> Processing: %s\n", cid))
    if(!file.exists(res_file)) { cat(sprintf("        WARN: Missing %s\n", res_file)); next }

    # LOAD & CHECK IDS
    res_df <- read.table(res_file, header=TRUE, sep="\t", quote="")
    if (grepl("^[0-9]+$", rownames(res_df)[1])) {
        id_col_idx <- which(colnames(res_df) %in% c("gene_id", "gene_name", "id", "GeneID"))
        if(length(id_col_idx) > 0) rownames(res_df) <- res_df[, id_col_idx[1]] else rownames(res_df) <- res_df[,1]
    }

    # Add symbols
    if(!"symbol" %in% colnames(res_df)) {
        clean_ids <- sub("\\..*", "", rownames(res_df))
        clean_ids <- gsub('"', '', clean_ids)
        syms <- mapIds(EnsDb.Hsapiens.v86, keys=clean_ids, column="SYMBOL", keytype="GENEID", multiVals="first")
        res_df$symbol <- ifelse(is.na(syms), clean_ids, syms)
    }

    # EXPORT RAW
    write.csv(head(res_df %>% arrange(padj), 50), paste0(gemini_dir, "/DE_Genes_", cid, ".csv"))

    # 1. VOLCANO (Size 9x8)
    p_vol <- EnhancedVolcano(res_df, lab=res_df$symbol, x='log2FoldChange', y='padj',
        title=paste0("Volcano: ", cid), subtitle=paste0(target, " vs ", ref),
        pCutoff=0.05, FCcutoff=1.0, pointSize=2.0, labSize=3.0,
        drawConnectors=TRUE, max.overlaps=30)
    save_plot(p_vol, paste0(output_prefix, "_Volcano_", cid, "_mqc"), 9, 8)

    # 2. GSEA
    res_df$rank <- res_df$log2FoldChange
    res_df <- res_df[!is.na(res_df$rank) & !is.infinite(res_df$rank),]
    gene_list <- sort(setNames(res_df$rank, res_df$symbol), decreasing=TRUE)

    gsea_out <- GSEA(gene_list, TERM2GENE=gmt, pvalueCutoff=1, verbose=FALSE, eps=1e-50)

    if(!is.null(gsea_out) && nrow(gsea_out) > 0) {
        if(nrow(gsea_out) > 50) gsea_out@result <- head(gsea_out@result[order(gsea_out@result$p.adjust), ], 50)
        write.csv(gsea_out@result[,c("ID", "NES", "p.adjust", "core_enrichment")], paste0(gemini_dir, "/GSEA_Results_", cid, ".csv"))

        top_up <- head(gsea_out@result[gsea_out@result$NES > 0, "ID"], 5)
        top_dn <- head(gsea_out@result[gsea_out@result$NES < 0, "ID"], 5)
        prompt_summary[[cid]] <- paste0("Contrast ", cid, ":\n- Top Up: ", paste(top_up, collapse=", "), "\n- Top Dn: ", paste(top_dn, collapse=", "), "\n")

        gsea_out <- pairwise_termsim(gsea_out)
        save_plot(dotplot(gsea_out, showCategory=15, split=".sign") + facet_grid(.~.sign), paste0(output_prefix, "_GSEA_Dot_", cid, "_mqc"), 12, 10)
        
        tryCatch({ save_plot(emapplot(gsea_out, showCategory=20, cex.params=list(category_label=0.7)), paste0(output_prefix, "_GSEA_Emap_", cid, "_mqc"), 14, 12) }, error=function(e) NULL)
    }
}

cat("LOG [5/5]: Generating Prompt...\n")
writeLines(paste0("Loaded Results from Pipeline:\n\n", paste(unlist(prompt_summary), collapse="\n")), paste0(dirname(output_prefix), "/GEMINI_PROMPTS.txt"))
