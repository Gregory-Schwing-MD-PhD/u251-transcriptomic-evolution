#!/usr/bin/env Rscript

# ==========================================
# MICROENVIRONMENT VISUALIZATION SUITE
# ==========================================
# This script generates publication-ready plots for Host (Rat) response.
# CRITICAL FEATURE: It accepts a "Clean Genes" list (Arg 6) to filter out
# background brain noise before generating heatmaps and pathway analysis.

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ComplexHeatmap)
    library(circlize)
    library(EnhancedVolcano)
    library(clusterProfiler)
    library(enrichplot)
    library(stringr)
})

# --- ARGUMENTS ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
    stop("Usage: Rscript plot_kitchen_sink_microenvironment.R <DESeq2_Res> <VST> <GMT> <OutPrefix> <Counts> <CleanGenesCSV>")
}

deseq_file   <- args[1]
vst_file     <- args[2]
gmt_file     <- args[3]
out_prefix   <- args[4]
counts_file  <- args[5]
clean_genes_file <- args[6]

# Create output directory if needed
out_dir <- dirname(out_prefix)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ==========================================
# 1. LOAD DATA & FILTERING
# ==========================================
cat("LOG: Loading Data...\n")
res_df <- fread(deseq_file)
vst_df <- fread(vst_file)
clean_genes <- fread(clean_genes_file)

# --- THE FILTER ---
# We restrict the universe of "significant" genes to ONLY those that passed subtraction
valid_genes <- clean_genes$Gene
cat(sprintf("LOG: Filtering analysis to %d validated LITT-response genes (Signal Subtracted).\n", length(valid_genes)))

# Mark genes in the main results table
# We add a column "Is_Validated" for the Volcano plot
res_df$Is_Validated <- ifelse(res_df$gene_id %in% valid_genes, "TRUE", "FALSE")

# Subset VST for Heatmaps (Only validated genes)
vst_mat <- as.matrix(vst_df[, -1, with=FALSE])
rownames(vst_mat) <- vst_df$gene_id
vst_clean <- vst_mat[rownames(vst_mat) %in% valid_genes, ]

# ==========================================
# 2. ENHANCED VOLCANO (With Validation Highlight)
# ==========================================
cat("LOG: Generating Volcano Plot...\n")

# We plot ALL genes, but we highlight only the validated ones in Red
# Non-validated high LFC genes (artifacts) will appear in Grey/Blue
keyvals <- ifelse(
    res_df$gene_id %in% valid_genes & res_df$log2FoldChange > 0, 'red',
    ifelse(res_df$gene_id %in% valid_genes & res_df$log2FoldChange < 0, 'blue',
    'grey'))

names(keyvals)[keyvals == 'red']  <- 'Validated Upregulated'
names(keyvals)[keyvals == 'blue'] <- 'Validated Downregulated'
names(keyvals)[keyvals == 'grey'] <- 'Background/Artifact'

pdf(paste0(out_prefix, "_Volcano_Validated.pdf"), width=10, height=8)
EnhancedVolcano(res_df,
    lab = res_df$gene_name,
    x = 'log2FoldChange',
    y = 'padj',
    title = 'LITT Response (Signal Subtracted)',
    subtitle = 'Red/Blue = Validated Candidates vs. Artifacts',
    pCutoff = 0.05,
    FCcutoff = 1.0,
    colCustom = keyvals,
    legendPosition = 'right',
    drawConnectors = TRUE,
    widthConnectors = 0.5
)
dev.off()

# ==========================================
# 3. HEATMAP (Validated Genes Only)
# ==========================================
cat("LOG: Generating Heatmap...\n")

if (nrow(vst_clean) > 2) {
    # Scale rows (Z-score)
    mat_scaled <- t(scale(t(vst_clean)))
    
    # Simple annotation (extract condition from column names if possible, else guess)
    # Assuming standard nf-core naming: sample_replicate
    col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

    pdf(paste0(out_prefix, "_Heatmap_CleanGenes.pdf"), width=10, height=length(valid_genes)/5 + 5)
    ht <- Heatmap(mat_scaled, 
        name = "Z-score",
        col = col_fun,
        show_row_names = length(valid_genes) < 100, # Only show names if list is short
        show_column_names = TRUE,
        cluster_columns = TRUE,
        column_title = "Validated Microenvironment Response"
    )
    draw(ht)
    dev.off()
} else {
    cat("WARNING: Not enough valid genes for heatmap.\n")
}

# ==========================================
# 4. PATHWAY ANALYSIS (ORA on Validated Genes)
# ==========================================
cat("LOG: Running Pathway Analysis (ORA)...\n")

# Load GMT
pathways <- read.gmt(gmt_file)

# For ORA, we test the Validated Genes against the background of ALL detected genes
universe <- res_df$gene_id
sig_genes_up <- clean_genes[clean_genes$LFC_LITT > 0]$Gene
sig_genes_down <- clean_genes[clean_genes$LFC_LITT < 0]$Gene

run_ora <- function(genes, title_suffix) {
    if (length(genes) > 0) {
        ego <- enricher(genes,
            TERM2GENE = pathways,
            universe = universe,
            pvalueCutoff = 0.1,
            qvalueCutoff = 0.2
        )
        
        if (!is.null(ego) && nrow(ego@result) > 0) {
            # Dotplot
            pdf(paste0(out_prefix, "_ORA_Dotplot_", title_suffix, ".pdf"), width=10, height=8)
            print(dotplot(ego, showCategory=20) + ggtitle(paste("Pathways:", title_suffix)))
            dev.off()
            
            # Save Table
            write.csv(ego@result, paste0(out_prefix, "_ORA_Table_", title_suffix, ".csv"))
        } else {
            cat(paste("LOG: No significant pathways for", title_suffix, "\n"))
        }
    }
}

run_ora(sig_genes_up, "Upregulated_Clean")
run_ora(sig_genes_down, "Downregulated_Clean")

cat("LOG: Analysis Complete.\n")
