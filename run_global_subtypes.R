#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# GLOBAL SUBTYPES & EVOLUTIONARY ANALYSIS v5.0
# ------------------------------------------------------------------------------
# - Focus: Trajectory (PCA Arrows), Plasticity (Entropy), High-Sensitivity Stats
# - Stats: Limma (eBayes), ANOVA, Bartlett, PERMANOVA
# - Output: Dynamic HTML Report + _mqc.png plots
# ------------------------------------------------------------------------------

set.seed(12345)

suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(ape); library(ggrepel)
    library(EnsDb.Hsapiens.v86); library(ComplexHeatmap); library(circlize)
    library(GSVA); library(tidyr); library(tibble)
    library(limma)
    library(vegan)
    library(car)
})

# ==============================================================================
# CONFIGURATION
# ==============================================================================
N_TOP_VARIABLE_GENES <- 500
PCA_POINT_SIZE <- 6
TEXT_BASE_SIZE <- 14

# ==============================================================================
# REPORTING ENGINE
# ==============================================================================
html_buffer <- character()

init_html <- function() {
    style <- "
    <style>
        body { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; color: #333; background: #fff; padding: 20px; }
        .report-container { max-width: 1200px; margin: 0 auto; }
        h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
        h2 { color: #2980b9; margin-top: 40px; border-left: 5px solid #2980b9; padding-left: 10px; }
        h3 { color: #16a085; margin-top: 30px; }
        table { width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 12px; }
        th { background: #ecf0f1; padding: 10px; text-align: left; border-bottom: 2px solid #bdc3c7; }
        td { padding: 8px; border-bottom: 1px solid #ecf0f1; }
        .sig { color: #c0392b; font-weight: bold; background-color: #f9e79f; }
        .img-box { text-align: center; margin: 30px 0; border: 1px solid #eee; padding: 10px; }
        .img-box img { max-width: 90%; height: auto; }
        .caption { margin-top: 10px; font-style: italic; color: #7f8c8d; }
        .citation-box { background: #f8f9fa; padding: 15px; margin-top: 50px; border-top: 2px solid #ddd; font-size: 12px; color: #555; }
    </style>
    <div class='report-container'>
    <h1>Evolutionary Subtype Analysis</h1>
    <p><strong>Overview:</strong> This report analyzes the trajectory of U251 evolution using PCA, Plasticity (Entropy) scoring, and High-Sensitivity Differential Analysis (Limma).</p>"
    html_buffer <<- c(html_buffer, style)
}

add_header <- function(txt) { html_buffer <<- c(html_buffer, paste0("<h2>", txt, "</h2>")) }

add_table <- function(df, title) {
    if(nrow(df) == 0) return()
    header <- paste0("<tr>", paste0("<th>", colnames(df), "</th>", collapse=""), "</tr>")
    rows <- apply(df, 1, function(r) {
        cells <- sapply(r, function(x) {
            val <- as.character(x)
            # Simple heuristic for highlighting significance strings or small numbers
            if(grepl("e-", val) || (suppressWarnings(!is.na(as.numeric(val))) && as.numeric(val) < 0.05)) {
                return(paste0("<td class='sig'>", val, "</td>"))
            }
            return(paste0("<td>", val, "</td>"))
        })
        paste0("<tr>", paste(cells, collapse=""), "</tr>")
    })
    html_buffer <<- c(html_buffer, paste0("<h3>", title, "</h3><table>", header, paste(rows, collapse=""), "</table>"))
}

add_image <- function(file_path, caption) {
    fname <- basename(file_path)
    html_buffer <<- c(html_buffer, paste0(
        "<div class='img-box'>",
        "<img src='", fname, "' alt='", caption, "'>",
        "<div class='caption'>", caption, "</div></div>"
    ))
}

add_citations <- function() {
    cites <- "
    <div class='citation-box'>
        <h3>References</h3>
        <ul>
            <li><strong>Verhaak et al. (2010):</strong> <em>Integrated genomic analysis identifies clinically relevant subtypes of glioblastoma...</em> Cancer Cell. Defines Classical, Mesenchymal, Proneural.</li>
            <li><strong>Neftel et al. (2019):</strong> <em>An Integrative Model of Cellular States...</em> Cell. Defines dynamic states (AC-like, MES-like, OPC-like, NPC-like).</li>
            <li><strong>Garofano et al. (2021):</strong> <em>Pathway-based classification of glioblastoma...</em> Nature Cancer. Defines metabolic subtypes (Mitochondrial, Glycolytic).</li>
            <li><strong>Limma (Smyth 2004):</strong> <em>Linear Models for Microarray Data.</em> Uses Empirical Bayes to stabilize variance for small sample sizes.</li>
        </ul>
    </div>"
    html_buffer <<- c(html_buffer, cites)
}

finish_html <- function(prefix) {
    html_buffer <<- c(html_buffer, "</div>")
    out_file <- paste0(dirname(prefix), "/Global_Subtype_Report.html")
    writeLines(html_buffer, out_file)
    cat(paste0("SUCCESS: Report written to ", out_file, "\n"))
}

# ==============================================================================
# UTILITY
# ==============================================================================
map_genes_to_symbols <- function(gene_ids, db = EnsDb.Hsapiens.v86) {
    clean_ids <- sub("\\..*", "", gene_ids)
    symbols <- mapIds(db, keys=clean_ids, column="SYMBOL", keytype="GENEID", multiVals="first")
    ifelse(is.na(symbols), gene_ids, symbols)
}

save_plot_mqc <- function(plot_obj, prefix, suffix, w=9, h=8) {
    base <- paste0(prefix, "_", suffix, "_mqc") # Enforce _mqc suffix
    tryCatch({
        pdf(paste0(base, ".pdf"), width=w, height=h)
        if (inherits(plot_obj, "Heatmap")) draw(plot_obj) else print(plot_obj)
        dev.off()
        png_out <- paste0(base, ".png")
        png(png_out, width=w, height=h, units="in", res=300)
        if (inherits(plot_obj, "Heatmap")) draw(plot_obj) else print(plot_obj)
        dev.off()
        return(png_out)
    }, error = function(e) return(NULL))
}

# ==============================================================================
# MAIN PIPELINE
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: script.R <results_dir> <out_prefix> <meta>")

results_dir   <- args[1]
output_prefix <- args[2]
meta_file     <- args[3]

# Ensure output directory exists
dir.create(dirname(output_prefix), showWarnings = FALSE, recursive = TRUE)
init_html()

cat("LOG [1/7]: Loading Data...\n")
vst_file <- file.path(results_dir, "tables/processed_abundance/all.vst.tsv")
mat_vst <- as.matrix(read.table(vst_file, header=TRUE, row.names=1, check.names=FALSE))
meta <- read.csv(meta_file, row.names=1)
common <- intersect(colnames(mat_vst), rownames(meta))
mat_vst <- mat_vst[, common]; meta <- meta[common, , drop=FALSE]

if("Classification" %in% colnames(meta)) {
    # Force order for trajectory logic
    meta$Classification <- factor(meta$Classification, levels=c("Culture_U2", "Primary_U2", "Recurrent_U2"))
}

mapped_syms <- map_genes_to_symbols(rownames(mat_vst))
mat_sym <- as.data.frame(mat_vst) %>% tibble::rownames_to_column("symbol") %>%
    mutate(symbol = mapped_syms) %>% dplyr::filter(!is.na(symbol)) %>% group_by(symbol) %>%
    summarise(across(everything(), mean)) %>% tibble::column_to_rownames("symbol") %>% as.matrix()

# ==============================================================================
# 2. GLOBAL STRUCTURE & TRAJECTORY (Original PCA Logic)
# ==============================================================================
add_header("1. Global Evolutionary Trajectory")
cat("LOG [2/7]: Calculating PCA Trajectory...\n")

top_var <- head(order(apply(mat_vst, 1, var), decreasing=TRUE), N_TOP_VARIABLE_GENES)
mat_sig <- mat_vst[top_var, ]
pca <- prcomp(t(mat_sig))
var_pc <- round(summary(pca)$importance[2, 1:2]*100, 1)

pcaData <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Sample=rownames(meta), Class=meta$Classification)

# --- CENTROID & ARROW LOGIC (Restored from Original Script) ---
centroids <- aggregate(cbind(PC1, PC2) ~ Class, data=pcaData, FUN=mean)
# Ensure order matches factor levels for correct arrow drawing
centroids <- centroids[match(levels(meta$Classification), centroids$Class), ]
centroids <- na.omit(centroids) # Safety

# Arrow data: Culture -> Primary -> Recurrent (based on factor order)
arrow_data <- data.frame(
    x_start = centroids$PC1[-nrow(centroids)],
    y_start = centroids$PC2[-nrow(centroids)],
    x_end = centroids$PC1[-1],
    y_end = centroids$PC2[-1]
)

# PERMANOVA Stats
perm_res <- tryCatch({ adonis2(dist(t(mat_sig)) ~ Classification, data = meta, permutations = 999) }, error=function(e) NULL)
sub_txt <- if(!is.null(perm_res)) paste0("PERMANOVA: p=", perm_res$`Pr(>F)`[1], ", R2=", round(perm_res$R2[1]*100,1), "%") else ""

grp_cols <- c("Culture_U2"="#1f77b4", "Primary_U2"="#ff7f0e", "Recurrent_U2"="#d62728")

p_pca <- ggplot(pcaData, aes(x=PC1, y=PC2)) +
    geom_segment(data=arrow_data, aes(x=x_start, y=y_start, xend=x_end, yend=y_end),
                 arrow=arrow(length=unit(0.4,"cm"), type="closed"), color="grey50", linewidth=1.5, inherit.aes=FALSE) +
    geom_point(aes(fill=Class, shape=Class), size=7, color="black", stroke=0.8) +
    ggrepel::geom_text_repel(aes(label=Sample), size=4, box.padding=0.5) +
    geom_label(data=centroids, aes(x=PC1, y=PC2, label=Class), fill="white", alpha=0.8, fontface="bold") +
    scale_fill_manual(values=grp_cols) +
    scale_shape_manual(values=c(21, 24, 22)) +
    labs(title="Evolutionary Trajectory (PCA)", subtitle=sub_txt, 
         x=paste0("PC1 (", var_pc[1], "%)"), y=paste0("PC2 (", var_pc[2], "%)")) +
    theme_bw(base_size=TEXT_BASE_SIZE) + theme(legend.position="bottom")

img_pca <- save_plot_mqc(p_pca, output_prefix, "PCA_Trajectory", w=10, h=9)
add_image(img_pca, "Figure 1: Evolutionary Trajectory. Arrows indicate the shift between centroids of Culture, Primary, and Recurrent samples.")

# ==============================================================================
# 3. SUBTYPE SCORING
# ==============================================================================
cat("LOG [3/7]: Scoring Subtypes...\n")
sigs <- list(
    "Verhaak_MES"=c("CHI3L1","CD44","VIM","RELB","STAT3"),
    "Verhaak_CL"=c("EGFR","AKT2","NOTCH3","CCND2"),
    "Verhaak_PN"=c("PDGFRA","IDH1","OLIG2","SOX2","TP53"),
    "Neftel_MES"=c("CHI3L1","CD44","ANXA1","VIM"),
    "Neftel_AC"=c("APOE","AQP4","CLU","S100B"),
    "Neftel_OPC"=c("OLIG1","OLIG2","PDGFRA","SOX10"),
    "Neftel_NPC"=c("DLL3","SOX4","TUBB3","DCX"),
    "Garofano_MTC"=c("SLC45A1","GOT2","IDH2","SDHA","CS"),
    "Garofano_GPM"=c("SLC2A1","HK2","PFKP","LDHA","GAPDH")
)
gsva_res <- gsva(mat_sym, sigs, method="gsva", kcdf="Gaussian", verbose=FALSE)
gsva_scaled <- t(scale(t(gsva_res)))

# ==============================================================================
# 4. PLASTICITY (ENTROPY) - "Are cells getting confused?"
# ==============================================================================
add_header("2. Plasticity (Entropy Analysis)")
cat("LOG [4/7]: Calculating Entropy...\n")

calc_entropy <- function(scores) {
    probs <- exp(scores) / sum(exp(scores))
    -sum(probs * log(probs))
}
meta$Plasticity <- apply(gsva_res, 2, calc_entropy)

# ANOVA for Plasticity
plast_aov <- summary(aov(Plasticity ~ Classification, data=meta))
plast_p <- plast_aov[[1]][["Pr(>F)"]][1]

p_plast <- ggplot(meta, aes(x=Classification, y=Plasticity, fill=Classification)) +
    geom_boxplot(alpha=0.6) + geom_jitter(width=0.1, size=3) +
    scale_fill_manual(values=grp_cols) +
    labs(title="Subtype Plasticity (Shannon Entropy)", 
         subtitle=paste0("Does state confusion change? ANOVA P = ", formatC(plast_p, digits=2))) +
    theme_bw(base_size=14)

img_plast <- save_plot_mqc(p_plast, output_prefix, "Plasticity_Entropy", w=7, h=6)
add_image(img_plast, "Figure 2: Plasticity Score. Higher entropy suggests cells are less differentiated (mixed subtype states).")

# ==============================================================================
# 5. HIGH-SENSITIVITY STATS (LIMMA & VARIANCE)
# ==============================================================================
add_header("3. Statistical Analysis (Limma & Variance)")
cat("LOG [5/7]: Running Limma & Variance Tests...\n")

# A. Limma (Trajectory)
design <- model.matrix(~0 + meta$Classification)
colnames(design) <- levels(meta$Classification)
cont.matrix <- makeContrasts(
    Cult_vs_Prim = Culture_U2 - Primary_U2,
    Rec_vs_Prim  = Recurrent_U2 - Primary_U2,
    Rec_vs_Cult  = Recurrent_U2 - Culture_U2,
    levels = design
)
fit <- lmFit(gsva_res, design)
fit <- contrasts.fit(fit, cont.matrix)
fit <- eBayes(fit)

limma_res <- data.frame()
for(coef in colnames(cont.matrix)) {
    tmp <- topTable(fit, coef=coef, number=Inf)
    tmp$Subtype <- rownames(tmp)
    tmp$Comparison <- coef
    limma_res <- rbind(limma_res, tmp)
}
# Table for anything with P < 0.15 (Trend or better)
trend_res <- limma_res %>% dplyr::filter(P.Value < 0.15) %>% 
    dplyr::select(Subtype, Comparison, logFC, P.Value, adj.P.Val) %>% arrange(P.Value)
add_table(trend_res, "Differential Trajectories (Limma P < 0.15)")

# B. Variance Tests (Levene/Bartlett) - Is heterogeneity lost?
var_df <- data.frame()
for(sig in rownames(gsva_res)) {
    scores <- gsva_res[sig, ]
    # Bartlett is good for normal data
    bp <- tryCatch(bartlett.test(scores ~ meta$Classification)$p.value, error=function(e) NA)
    var_df <- rbind(var_df, data.frame(Subtype=sig, Bartlett_Variance_P=formatC(bp, format="e", digits=2)))
}
add_table(var_df, "Heterogeneity Tests (Variance Differences)")

# ==============================================================================
# 6. GLOBAL HEATMAP & RECOVERY
# ==============================================================================
add_header("4. Subtype Landscape")
cat("LOG [6/7]: Generating Heatmap...\n")

# Annotate heatmap with Plasticity
ha <- HeatmapAnnotation(
    Plasticity = anno_barplot(meta$Plasticity, height = unit(2, "cm")),
    Group = meta$Classification,
    col = list(Group = grp_cols)
)

ht <- Heatmap(gsva_scaled, name="Z-Score", 
    top_annotation = ha,
    column_split = meta$Classification,
    cluster_rows = TRUE, cluster_columns = TRUE,
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    column_title = "Global Subtype Landscape")

img_ht <- save_plot_mqc(ht, output_prefix, "Subtype_Heatmap", h=10)
add_image(img_ht, "Figure 3: Global Subtype Heatmap. Top barplot shows plasticity (entropy) per sample.")

# ==============================================================================
# 7. FINALIZE
# ==============================================================================
add_citations()
finish_html(output_prefix)
cat("\n=== PIPELINE COMPLETE ===\n")
