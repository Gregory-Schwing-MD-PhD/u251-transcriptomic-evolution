#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# EVOLUTIONARY KITCHEN SINK - GLOBAL ANALYSIS & HTML REPORTING
# ------------------------------------------------------------------------------
# - Focus: Global structure (PCA/MDS) & Comprehensive Subtype Analysis
# - Stats: PERMANOVA (Global) & Kruskal-Wallis (Subtype-specific)
# - Output: PDF/PNG Plots + Dynamic HTML Report (Mimicking v4 engine)
# ------------------------------------------------------------------------------

set.seed(12345)

suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(ape); library(ggrepel)
    library(EnsDb.Hsapiens.v86); library(ComplexHeatmap); library(circlize)
    library(GSVA); library(tidyr); library(tibble)
    if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan", repos="https://cloud.r-project.org")
    library(vegan)
    library(knitr) # For kable (optional, but good for text tables)
})

# ==============================================================================
# CONFIGURATION
# ==============================================================================
N_TOP_VARIABLE_GENES <- 500
PCA_POINT_SIZE <- 6
TEXT_BASE_SIZE <- 14

# ==============================================================================
# REPORTING ENGINE (Ported & Adapted from run_pathways_drugs_v4.R)
# ==============================================================================
html_buffer <- character()

init_html <- function() {
    style <- "
    <style>
        body { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; color: #333; background-color: #f4f7f6; margin: 0; padding: 20px; }
        .mqc-custom-content-section { max-width: 1200px; margin: 0 auto; background: #fff; padding: 30px; box-shadow: 0 0 15px rgba(0,0,0,0.1); border-radius: 5px; }
        h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 15px; }
        h2 { color: #34495e; margin-top: 30px; border-left: 5px solid #3498db; padding-left: 10px; }
        .stat-block { background: #e8f6f3; border: 1px solid #d4efdf; padding: 15px; border-radius: 4px; margin-bottom: 20px; }
        .stat-value { font-weight: bold; color: #16a085; font-size: 1.1em; }
        
        /* Table Styles */
        .subtype-table { width: 100%; border-collapse: collapse; margin-top: 15px; font-size: 13px; }
        .subtype-table th { background-color: #34495e; color: #fff; text-align: left; padding: 10px; }
        .subtype-table td { border-bottom: 1px solid #eee; padding: 8px; }
        .subtype-table tr:hover { background-color: #f9f9f9; }
        
        /* Significance Highlights */
        .sig-green { color: #27ae60; font-weight: bold; background-color: #eafaf1; padding: 3px 6px; border-radius: 3px; }
        .sig-grey { color: #95a5a6; font-style: italic; }
        
        /* Reference Section */
        .ref-box { background: #f8f9fa; border: 1px solid #eee; padding: 15px; margin-top: 30px; border-radius: 5px; font-size: 13px; }
        .ref-title { font-weight: bold; color: #2980b9; margin-bottom: 5px; display: block; }
        .ref-link { color: #3498db; text-decoration: none; }
        
        /* Image Container */
        .img-container { text-align: center; margin: 20px 0; padding: 10px; border: 1px solid #eee; background: #fafafa; }
        .img-container img { max-width: 100%; height: auto; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }
        .img-caption { margin-top: 10px; color: #7f8c8d; font-style: italic; font-size: 0.9em; }
    </style>
    <div class='mqc-custom-content-section'>
    <h1>Global Transcriptomic Analysis Report</h1>
    <p><strong>Generated:</strong> "
    
    html_buffer <<- c(html_buffer, paste0(style, Sys.Date(), "</p>"))
}

add_stat_block <- function(title, p_val, r2=NULL) {
    p_class <- ifelse(p_val < 0.05, "sig-green", "sig-grey")
    r2_html <- if(!is.null(r2)) paste0(" | Variance Explained: <strong>", r2, "%</strong>") else ""
    
    html <- paste0(
        "<div class='stat-block'>",
        "<strong>", title, "</strong><br>",
        "P-value: <span class='", p_class, "'>", formatC(p_val, format="e", digits=2), "</span>",
        r2_html,
        "</div>"
    )
    html_buffer <<- c(html_buffer, html)
}

add_subtype_table <- function(df) {
    rows <- apply(df, 1, function(r) {
        p_val <- as.numeric(r['P_Value'])
        sig_class <- ifelse(p_val < 0.05, "sig-green", "sig-grey")
        sig_label <- ifelse(p_val < 0.05, "Significant", "NS")
        
        paste0(
            "<tr>",
            "<td><strong>", r['Family'], "</strong></td>",
            "<td>", r['Subtype'], "</td>",
            "<td><span class='", sig_class, "'>", formatC(p_val, format="e", digits=2), "</span></td>",
            "<td>", sig_label, "</td>",
            "</tr>"
        )
    })
    
    table_html <- paste0(
        "<table class='subtype-table'>",
        "<tr><th>Family</th><th>Subtype Signature</th><th>Kruskal-Wallis P-value</th><th>Status</th></tr>",
        paste(rows, collapse=""),
        "</table>"
    )
    html_buffer <<- c(html_buffer, table_html)
}

add_image <- function(filepath, caption) {
    rel_path <- basename(filepath)
    html <- paste0(
        "<div class='img-container'>",
        "<img src='", rel_path, "' alt='", caption, "'>",
        "<div class='img-caption'>", caption, "</div>",
        "</div>"
    )
    html_buffer <<- c(html_buffer, html)
}

add_references <- function() {
    refs <- "
    <div class='ref-box'>
        <h3>Subtype Definitions & References</h3>
        
        <span class='ref-title'>1. Verhaak et al. (The Classics)</span>
        Paper: <a href='https://doi.org/10.1016/j.ccr.2009.12.020' class='ref-link' target='_blank'>Cancer Cell, 2010</a><br>
        <em>Includes: Classical (EGFR), Mesenchymal (NF1/Immune), Proneural (PDGFRA/IDH1).</em><br><br>
        
        <span class='ref-title'>2. Neftel et al. (Cellular States)</span>
        Paper: <a href='https://doi.org/10.1016/j.cell.2019.06.024' class='ref-link' target='_blank'>Cell, 2019</a><br>
        <em>Defines dynamic states: AC-like (Astrocytic), MES-like (Mesenchymal), NPC-like, OPC-like.</em><br><br>
        
        <span class='ref-title'>3. Garofano et al. (Metabolic Subtypes)</span>
        Paper: <a href='https://doi.org/10.1038/s43018-020-00159-4' class='ref-link' target='_blank'>Nature Cancer, 2021</a><br>
        <em>Focuses on metabolism: Mitochondrial (MTC), Glycolytic (GPM), Neuronal (NEU), Proliferative (PPR).</em>
    </div>
    "
    html_buffer <<- c(html_buffer, refs)
}

finish_html <- function(output_prefix) {
    html_buffer <<- c(html_buffer, "</div>")
    out_file <- paste0(dirname(output_prefix), "/Global_Subtype_Report_mqc.html")
    writeLines(html_buffer, out_file)
    cat(paste0("SUCCESS: Report generated at ", out_file, "\n"))
}

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================
map_genes_to_symbols <- function(gene_ids, db = EnsDb.Hsapiens.v86) {
    clean_ids <- sub("\\..*", "", gene_ids)
    symbols <- mapIds(db, keys=clean_ids, column="SYMBOL", keytype="GENEID", multiVals="first")
    ifelse(is.na(symbols), gene_ids, symbols)
}

save_plot <- function(plot_obj, filename_base, w=9, h=8) {
    # Saves PDF for publication, PNG for HTML report
    tryCatch({
        pdf(paste0(filename_base, ".pdf"), width=w, height=h)
        if (inherits(plot_obj, "Heatmap")) draw(plot_obj) else print(plot_obj)
        dev.off()
        
        png_file <- paste0(filename_base, ".png")
        png(png_file, width=w, height=h, units="in", res=300)
        if (inherits(plot_obj, "Heatmap")) draw(plot_obj) else print(plot_obj)
        dev.off()
        
        return(png_file) # Return filename for HTML embedding
    }, error = function(e) {
        cat(paste0("ERROR: Failed to save ", basename(filename_base), ": ", e$message, "\n"))
        return(NULL)
    })
}

# ==============================================================================
# MAIN PIPELINE
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: script.R <results_dir> <out_prefix> <meta>")

results_dir   <- args[1]
output_prefix <- args[2]
meta_file     <- args[3]

# Initialize Report
init_html()

cat("LOG [1/6]: Loading Data...\n")
vst_file <- file.path(results_dir, "tables/processed_abundance/all.vst.tsv")
mat_vst <- as.matrix(read.table(vst_file, header=TRUE, row.names=1, check.names=FALSE))
meta <- read.csv(meta_file, row.names=1)
common <- intersect(colnames(mat_vst), rownames(meta))
mat_vst <- mat_vst[, common]; meta <- meta[common, , drop=FALSE]

if("Classification" %in% colnames(meta)) {
    meta$Classification <- factor(meta$Classification, levels=c("Culture_U2", "Primary_U2", "Recurrent_U2"))
}

# Prepare symbols for GSVA
mapped_syms <- map_genes_to_symbols(rownames(mat_vst))
mat_sym <- as.data.frame(mat_vst) %>% tibble::rownames_to_column("symbol") %>%
    mutate(symbol = mapped_syms) %>% filter(!is.na(symbol)) %>% group_by(symbol) %>%
    summarise(across(everything(), mean)) %>% tibble::column_to_rownames("symbol") %>% as.matrix()

# ==============================================================================
# DEFINE SIGNATURES
# ==============================================================================
verhaak_sigs <- list(
    "Verhaak_Mesenchymal" = c("CHI3L1", "CD44", "VIM", "RELB", "STAT3", "MET", "TRADD", "MMP9", "TIMP1"),
    "Verhaak_Classical"   = c("EGFR", "AKT2", "NOTCH3", "JAG1", "CCND2", "F3", "PDGFA", "NES"),
    "Verhaak_Proneural"   = c("PDGFRA", "IDH1", "OLIG2", "SOX2", "NKX2-2", "OLIG1", "TP53"),
    "Verhaak_Neural"      = c("NEFL", "GABRA1", "SYT1", "SLC12A5", "MBP", "GABRG2")
)
garofano_sigs <- list(
    "Garofano_Mitochondrial" = c("SLC45A1", "GOT2", "ACO2", "IDH2", "IDH3A", "OGDH", "SUCLA2", "SDHA", "SDHB", "SDHC", "SDHD", "FH", "MDH2", "CS", "DLAT"),
    "Garofano_Glycolytic"    = c("SLC2A1", "HK2", "PFKP", "ALDOA", "GAPDH", "PGK1", "ENO1", "PKM", "LDHA", "SLC16A3", "PDK1"),
    "Garofano_Neuronal"      = c("SYP", "DLG4", "GRIN1", "GRIN2B", "GABRA1", "NEFL", "STMN2", "TUBB3", "MAP2"),
    "Garofano_Proliferative" = c("MKI67", "TOP2A", "PCNA", "CCNB1", "CCNB2", "CDK1", "AURKA", "AURKB", "BIRC5")
)
neftel_sigs <- list(
    "Neftel_AC_Like"  = c("APOE", "AQP4", "CLU", "CST3", "GPX3", "S100B", "SLC1A2", "SLC1A3"),
    "Neftel_MES_Like" = c("CHI3L1", "CD44", "ANXA1", "ANXA2", "C3", "CD74", "TIMP1", "VIM"),
    "Neftel_NPC_Like" = c("DLL3", "SOX4", "SOX11", "STMN1", "TAGLN3", "TUBB3", "DCX"),
    "Neftel_OPC_Like" = c("OLIG1", "OLIG2", "PDGFRA", "SOX10", "NKX2-2", "CSPG4", "PTPRZ1")
)
all_sigs <- c(verhaak_sigs, garofano_sigs, neftel_sigs)

# ==============================================================================
# GLOBAL STATS & PLOTS
# ==============================================================================
cat("LOG [2/6]: Global Stats & Plots...\n")
html_buffer <<- c(html_buffer, "<h2>1. Global Sample Structure</h2>")

top_var <- head(order(apply(mat_vst, 1, var), decreasing=TRUE), N_TOP_VARIABLE_GENES)
mat_sig <- mat_vst[top_var, ]
dist_mat <- dist(t(mat_sig))

# PERMANOVA
perm_res <- tryCatch({ adonis2(dist_mat ~ Classification, data = meta, permutations = 999) }, error = function(e) NULL)
if(!is.null(perm_res)) {
    perm_p <- perm_res$`Pr(>F)`[1]
    perm_r2 <- round(perm_res$R2[1] * 100, 1)
    add_stat_block("Global PERMANOVA (Separation by Group)", perm_p, perm_r2)
    stat_subtitle <- paste0("PERMANOVA: p=", perm_p, ", R2=", perm_r2, "%")
} else {
    stat_subtitle <- "PERMANOVA: Failed"
}

# Tree
phylo_tree <- as.phylo(hclust(dist_mat, method="ward.D2"))
if("C2B" %in% phylo_tree$tip.label) phylo_tree <- root(phylo_tree, outgroup="C2B", resolve.root=TRUE)
group_colors <- c("Culture" = "#1f77b4", "Primary" = "#ff7f0e", "Recurrent" = "#d62728")
tip_colors <- if("Classification" %in% colnames(meta)) group_colors[as.character(meta[phylo_tree$tip.label, "Classification"])] else "black"
tip_colors[is.na(tip_colors)] <- "black"

pdf(paste0(output_prefix, "_Phylogenetic_Tree.pdf"), width=10, height=8)
plot(phylo_tree, main="Phylogenetic Tree", sub=stat_subtitle, type="phylogram", edge.width=2, tip.color=tip_colors, label.offset=0.5, cex=1.2)
dev.off()
# Re-save as PNG for HTML (since base plot is tricky to capture otherwise)
png_tree <- paste0(output_prefix, "_Phylogenetic_Tree.png")
png(png_tree, width=10, height=8, units="in", res=300)
plot(phylo_tree, main="Phylogenetic Tree", sub=stat_subtitle, type="phylogram", edge.width=2, tip.color=tip_colors, label.offset=0.5, cex=1.2)
dev.off()
add_image(png_tree, "Figure 1: Phylogenetic Tree (Ward's D2 Clustering)")

# PCA
pca <- prcomp(t(mat_sig))
pcaData <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Sample=rownames(meta), Class=meta$Classification)
var_expl <- round(summary(pca)$importance[2, 1:2] * 100, 1)

p_pca <- ggplot(pcaData, aes(x=PC1, y=PC2, fill=Class, shape=Class)) +
    geom_point(size=PCA_POINT_SIZE, color="black", stroke=0.8) +
    ggrepel::geom_text_repel(aes(label=Sample), box.padding=0.5, size=4) +
    scale_shape_manual(values=c(21, 24, 22)) + scale_fill_manual(values=group_colors) +
    labs(title="PCA Global", subtitle=stat_subtitle, x=paste0("PC1 (", var_expl[1], "%)"), y=paste0("PC2 (", var_expl[2], "%)")) +
    theme_bw(base_size=TEXT_BASE_SIZE)
png_pca <- save_plot(p_pca, paste0(output_prefix, "_PCA_Global"))
add_image(png_pca, "Figure 2: Principal Component Analysis (PCA)")

# ==============================================================================
# SUBTYPE ANALYSIS
# ==============================================================================
cat("LOG [3/6]: Subtype Analysis & Stats...\n")
html_buffer <<- c(html_buffer, "<h2>2. GBM Subtype & State Scores</h2>")

gsva_res <- gsva(mat_sym, all_sigs, method="gsva", kcdf="Gaussian", verbose=FALSE)
gsva_norm <- apply(gsva_res, 2, function(x) round((x-min(x))/(max(x)-min(x))*100, 0))

# Kruskal-Wallis & Table Generation
stats_df <- data.frame(Family=character(), Subtype=character(), P_Value=numeric(), stringsAsFactors=FALSE)

subtype_pvals <- apply(gsva_res, 1, function(scores) {
    tryCatch({ kruskal.test(scores ~ meta$Classification)$p.value }, error = function(e) NA)
})

# Populate Table Data
for(nm in names(subtype_pvals)) {
    parts <- strsplit(nm, "_")[[1]]
    fam <- parts[1]
    sub <- paste(parts[-1], collapse=" ")
    stats_df <- rbind(stats_df, data.frame(Family=fam, Subtype=sub, P_Value=subtype_pvals[nm]))
}
add_subtype_table(stats_df)

# Heatmap
new_rownames <- paste0(names(subtype_pvals), " (p=", format(subtype_pvals, digits=2), ")")
rownames(gsva_norm) <- new_rownames
row_split_vec <- c(rep("Verhaak", 4), rep("Garofano", 4), rep("Neftel", 4))

ht_sub <- Heatmap(gsva_norm, name = "Score", column_split = meta$Classification,
    row_split = row_split_vec, cluster_rows = FALSE, cluster_columns = TRUE,
    col = colorRamp2(c(0, 50, 100), c("white", "yellow", "red")),
    row_names_gp = gpar(fontsize = 10), column_title = "GBM Subtype Scores")
png_ht <- save_plot(ht_sub, paste0(output_prefix, "_Subtype_Heatmap"), h=12)
add_image(png_ht, "Figure 3: Subtype Signature Heatmap (Normalized Scores)")

# ==============================================================================
# FINALIZE
# ==============================================================================
add_references()
finish_html(output_prefix)
cat("\n=== PIPELINE COMPLETE ===\n")
