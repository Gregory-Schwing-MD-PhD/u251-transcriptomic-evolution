#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# THE KITCHEN SINK: MASTER VISUALIZATION SUITE (v5 - Human)
# ------------------------------------------------------------------------------
# Features:
# 1. Robust ID Mapping (Ensembl -> Symbol) using org.Hs.eg.db
# 2. Enhanced Volcano Plots (labeled with Symbols)
# 3. GSEA & GSVA (Pathway Analysis)

suppressPackageStartupMessages({
    library(DESeq2)
    library(clusterProfiler)
    library(enrichplot)
    library(GSVA)
    library(GSEABase)
    library(ggplot2)
    library(ggnewscale)
    library(EnhancedVolcano)
    library(ComplexHeatmap)
    library(circlize)
    library(stringr)
    library(data.table)
    library(org.Hs.eg.db) # Human Database
})

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript plot_kitchen_sink.R <deseq_res_file> <vst_counts_file> <gmt_file> <output_prefix>")
}

deseq_file  <- args[1]
vst_file    <- args[2]
gmt_file    <- args[3]
out_prefix  <- args[4]

# Create output directory
dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. ID MAPPING FUNCTION (HUMAN)
# ==============================================================================
map_ids_to_symbols <- function(ids) {
    # Remove version numbers (ENSG0001.2 -> ENSG0001)
    clean_ids <- sub("\\..*", "", ids)
    
    # Map using the Human Database
    syms <- mapIds(org.Hs.eg.db, keys = clean_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
    
    # Fill NAs with original ID (fallback)
    syms[is.na(syms)] <- ids[is.na(syms)]
    
    return(make.unique(as.character(syms)))
}

# ==============================================================================
# 2. LOAD & MAP DATA
# ==============================================================================
cat("LOG: Loading and Mapping Data...\n")

# A. DESeq2 Results
res <- data.table::fread(deseq_file)
# Detect ID column (usually gene_id or the first column)
if ("gene_id" %in% colnames(res)) {
    res$symbol <- map_ids_to_symbols(res$gene_id)
} else {
    # If no gene_id col, assume rownames are IDs (but fread doesn't preserve rownames well)
    # Ideally, input has gene_id. If not, we take the first column.
    res$symbol <- map_ids_to_symbols(res[[1]])
}

# B. VST Counts
vst_dt <- data.table::fread(vst_file)
colnames(vst_dt)[1] <- "gene_id"
vst_dt$symbol <- map_ids_to_symbols(vst_dt$gene_id)

# Aggregate VST by Symbol (Mean) - Handle duplicates if any
numeric_cols <- setdiff(colnames(vst_dt), c("gene_id", "symbol"))
vst_aggr <- vst_dt[, lapply(.SD, mean), by=symbol, .SDcols=numeric_cols]
vst_mat <- as.matrix(vst_aggr[, -1, with=FALSE])
rownames(vst_mat) <- vst_aggr$symbol

# ==============================================================================
# 3. VISUALIZATIONS
# ==============================================================================

# --- Volcano Plot ---
cat("LOG: Generating Volcano Plot...\n")
p_vol <- EnhancedVolcano(
    res,
    lab = res$symbol,
    x = 'log2FoldChange',
    y = 'pvalue',
    title = 'Differential Expression',
    subtitle = 'Recurrent vs Primary',
    pCutoff = 0.05,
    FCcutoff = 1.0,
    pointSize = 2.0,
    labSize = 4.0,
    legendPosition = 'right',
    max.overlaps = 20,
    drawConnectors = TRUE,
    widthConnectors = 0.5
)
ggsave(paste0(out_prefix, "_Volcano_mqc.pdf"), p_vol, width=10, height=8)

# --- GSEA ---
cat("LOG: Running GSEA...\n")
# Rank metric: Signed -log10(pvalue)
res$rank_metric <- sign(res$log2FoldChange) * -log10(res$pvalue)
res <- res[!is.na(rank_metric) & !is.infinite(rank_metric)]
res <- res[order(abs(rank_metric), decreasing = TRUE)]

# Remove duplicates (keep highest rank)
res <- res[!duplicated(symbol)]
gene_list <- setNames(res$rank_metric, res$symbol)
gene_list <- sort(gene_list, decreasing = TRUE)

gmt <- read.gmt(gmt_file)
gsea_res <- GSEA(gene_list, TERM2GENE = gmt, pvalueCutoff = 1.0, verbose = FALSE, eps=1e-50)

if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
    top_gsea <- gsea_res
    # Top 50 by significance
    if (nrow(top_gsea) > 50) top_gsea@result <- head(top_gsea@result[order(top_gsea@result$p.adjust), ], 50)
    top_gsea <- pairwise_termsim(top_gsea)

    # Dotplot
    p_dot <- dotplot(top_gsea, showCategory=15, split=".sign") + facet_grid(.~.sign)
    ggsave(paste0(out_prefix, "_GSEA_Dotplot_mqc.pdf"), p_dot, width=12, height=8)
    
    # Cnet Plot
    p_cnet <- cnetplot(top_gsea, categorySize="pvalue", foldChange=gene_list, circular=TRUE, color.params=list(edge=TRUE))
    ggsave(paste0(out_prefix, "_GSEA_Network_mqc.pdf"), p_cnet, width=12, height=12)
}

# --- GSVA Heatmap ---
cat("LOG: Running GSVA...\n")
gene_sets_list <- split(gmt$gene, gmt$term)
gsva_res <- gsva(vst_mat, gene_sets_list, method="gsva", kcdf="Gaussian", verbose=FALSE)

pathway_var <- apply(gsva_res, 1, var)
top_pathways <- names(sort(pathway_var, decreasing=TRUE))[1:30] # Top 30 variable pathways

ht <- Heatmap(gsva_res[top_pathways, ], 
              name="GSVA", 
              column_title="Pathway Activity", 
              col=colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")))

pdf(paste0(out_prefix, "_GSVA_Heatmap_mqc.pdf"), width=10, height=8)
draw(ht)
dev.off()

cat("LOG: Analysis Complete.\n")
