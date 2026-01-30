#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# GLOBAL SUBTYPES & EVOLUTIONARY ANALYSIS v7.0 (Self-Contained)
# ------------------------------------------------------------------------------
# - Focus: Trajectory, Plasticity, High-Sensitivity Stats
# - Fixes: Base64 Image Embedding, Raw Data Tables, Tree Generation
# ------------------------------------------------------------------------------

set.seed(12345)

suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(ape); library(ggrepel)
    library(EnsDb.Hsapiens.v86); library(ComplexHeatmap); library(circlize)
    library(GSVA); library(tidyr); library(tibble)
    library(limma)
    library(vegan)
    library(car)
    library(base64enc) # Critical for self-contained reports
})

# ==============================================================================
# CONFIGURATION
# ==============================================================================
N_TOP_VARIABLE_GENES <- 500
PCA_POINT_SIZE <- 6
TEXT_BASE_SIZE <- 14

# ==============================================================================
# REPORTING ENGINE (Base64 Embedded)
# ==============================================================================
html_buffer <- character()

init_html <- function() {
    style <- "
    <style>
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; color: #333; background: #fff; padding: 40px; max-width: 1200px; margin: 0 auto; }
        h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 15px; }
        h2 { color: #2980b9; margin-top: 50px; border-left: 6px solid #2980b9; padding-left: 15px; background: #f4f6f7; padding: 10px; }
        h3 { color: #16a085; margin-top: 30px; border-bottom: 1px solid #eee; padding-bottom: 5px;}
        p { line-height: 1.6; color: #555; margin-bottom: 15px; }
        
        /* Tables */
        table { width: 100%; border-collapse: collapse; margin-top: 15px; font-size: 12px; border: 1px solid #ddd; }
        th { background: #34495e; color: #fff; padding: 12px; text-align: left; }
        td { padding: 8px; border-bottom: 1px solid #eee; }
        tr:nth-child(even) { background: #f9f9f9; }
        .sig { color: #c0392b; font-weight: bold; background-color: #fdebd0; }
        
        /* Images */
        .img-box { text-align: center; margin: 40px 0; border: 1px solid #ddd; padding: 20px; background: #fafafa; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }
        .img-box img { max-width: 95%; height: auto; }
        .caption { margin-top: 15px; font-style: italic; color: #7f8c8d; font-size: 1em; font-weight: 500; }
        
        /* Misc */
        .citation-box { background: #eaf2f8; padding: 20px; margin-top: 60px; border-top: 4px solid #3498db; font-size: 13px; color: #444; border-radius: 5px; }
        .stat-explainer { background: #fff8e1; border-left: 4px solid #f1c40f; padding: 15px; font-size: 0.95em; margin-bottom: 20px; color: #7f6000; }
        .llm-box { background: #2c3e50; color: #ecf0f1; border-radius: 5px; padding: 20px; margin-top: 40px; font-family: monospace; font-size: 12px; white-space: pre-wrap; }
        .scroll-table { overflow-x: auto; display: block; max-height: 500px; }
    </style>
    <h1>Global Evolutionary Transcriptomic Report</h1>
    <p><strong>Generated:</strong> "
    html_buffer <<- c(html_buffer, paste0(style, Sys.Date(), "</p>"))
}

add_header <- function(txt) { html_buffer <<- c(html_buffer, paste0("<h2>", txt, "</h2>")) }

add_explainer <- function(title, text) {
    html_buffer <<- c(html_buffer, paste0("<div class='stat-explainer'><strong>", title, ":</strong> ", text, "</div>"))
}

add_table <- function(df, title, scroll=FALSE) {
    if(nrow(df) == 0) {
        html_buffer <<- c(html_buffer, paste0("<h3>", title, "</h3><p>No results found.</p>"))
        return()
    }
    
    # Format numeric columns
    df_fmt <- df
    for(i in 1:ncol(df)) {
        if(is.numeric(df[,i])) df_fmt[,i] <- formatC(df[,i], format="f", digits=3)
    }
    
    header <- paste0("<tr>", paste0("<th>", colnames(df), "</th>", collapse=""), "</tr>")
    rows <- apply(df_fmt, 1, function(r) {
        cells <- sapply(r, function(x) {
            val <- as.character(x)
            # Heuristic for highlighting
            if(grepl("e-", val) || (suppressWarnings(!is.na(as.numeric(val))) && as.numeric(val) < 0.05 && as.numeric(val) > 0)) {
                return(paste0("<td class='sig'>", val, "</td>"))
            }
            return(paste0("<td>", val, "</td>"))
        })
        paste0("<tr>", paste(cells, collapse=""), "</tr>")
    })
    
    tbl_class <- if(scroll) "scroll-table" else ""
    html_buffer <<- c(html_buffer, paste0("<h3>", title, "</h3><div class='", tbl_class, "'><table>", header, paste(rows, collapse=""), "</table></div>"))
}

add_image_base64 <- function(file_path, caption) {
    if(!file.exists(file_path)) {
        html_buffer <<- c(html_buffer, paste0("<p style='color:red'>Error: Image not found - ", basename(file_path), "</p>"))
        return()
    }
    # Encode image to Base64
    b64 <- base64enc::dataURI(file = file_path, mime = "image/png")
    html_buffer <<- c(html_buffer, paste0(
        "<div class='img-box'>",
        "<img src='", b64, "' alt='", caption, "'>",
        "<div class='caption'>", caption, "</div></div>"
    ))
}

add_llm_summary <- function(txt) {
    html_buffer <<- c(html_buffer, paste0("<h3>AI Summary Context</h3><div class='llm-box'><strong>Prompt Context:</strong><br>", txt, "</div>"))
}

add_citations <- function() {
    cites <- "
    <div class='citation-box'>
        <h3>References & Methodology</h3>
        <ul>
            <li><strong>Verhaak et al. (2010):</strong> <em>Cancer Cell.</em> (Classical, Mesenchymal, Proneural).</li>
            <li><strong>Neftel et al. (2019):</strong> <em>Cell.</em> (AC-like, MES-like, OPC-like, NPC-like).</li>
            <li><strong>Garofano et al. (2021):</strong> <em>Nature Cancer.</em> (Mitochondrial, Glycolytic).</li>
            <li><strong>Limma:</strong> Linear Models for Microarray Data (Smyth 2004). Empirical Bayes moderation.</li>
            <li><strong>ANOVA & Bartlett:</strong> Parametric tests for mean and variance differences.</li>
        </ul>
    </div>"
    html_buffer <<- c(html_buffer, cites)
}

finish_html <- function(prefix) {
    out_file <- paste0(dirname(prefix), "/U251_Global_Subtypes_Report_Standalone.html")
    writeLines(html_buffer, out_file)
    cat(paste0("SUCCESS: Self-contained report written to ", out_file, "\n"))
}

# ==============================================================================
# UTILITY
# ==============================================================================
map_genes_to_symbols <- function(gene_ids, db = EnsDb.Hsapiens.v86) {
    clean_ids <- sub("\\..*", "", gene_ids)
    symbols <- mapIds(db, keys=clean_ids, column="SYMBOL", keytype="GENEID", multiVals="first")
    ifelse(is.na(symbols), gene_ids, symbols)
}

safe_format <- function(x) {
    if (is.null(x) || length(x) == 0) return("NA")
    if (is.na(x)) return("NA")
    if (!is.numeric(x)) return(as.character(x))
    return(formatC(x, format="e", digits=2))
}

# ==============================================================================
# MAIN PIPELINE
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: script.R <results_dir> <out_prefix> <meta>")

results_dir   <- args[1]
output_prefix <- args[2]
meta_file     <- args[3]

dir.create(dirname(output_prefix), showWarnings = FALSE, recursive = TRUE)
init_html()

cat("LOG [1/8]: Loading Data...\n")
vst_file <- file.path(results_dir, "tables/processed_abundance/all.vst.tsv")
mat_vst <- as.matrix(read.table(vst_file, header=TRUE, row.names=1, check.names=FALSE))
meta <- read.csv(meta_file, row.names=1)
common <- intersect(colnames(mat_vst), rownames(meta))
mat_vst <- mat_vst[, common]; meta <- meta[common, , drop=FALSE]

if("Classification" %in% colnames(meta)) {
    meta$Classification <- factor(meta$Classification, levels=c("Culture_U2", "Primary_U2", "Recurrent_U2"))
}

mapped_syms <- map_genes_to_symbols(rownames(mat_vst))
mat_sym <- as.data.frame(mat_vst) %>% tibble::rownames_to_column("symbol") %>%
    mutate(symbol = mapped_syms) %>% dplyr::filter(!is.na(symbol)) %>% group_by(symbol) %>%
    summarise(across(everything(), mean)) %>% tibble::column_to_rownames("symbol") %>% as.matrix()

# ==============================================================================
# 2. GLOBAL STRUCTURE (TREE & PCA)
# ==============================================================================
add_header("1. Global Evolutionary Trajectory")
add_explainer("Phylogenetic Tree", "Visualizes the distance between samples based on whole-transcriptome expression. Branches closer together are more similar.")

cat("LOG [2/8]: Calculating PCA & Tree...\n")

top_var <- head(order(apply(mat_vst, 1, var), decreasing=TRUE), N_TOP_VARIABLE_GENES)
mat_sig <- mat_vst[top_var, ]

# --- PHYLOGENETIC TREE (Explicit Base Save) ---
phylo_tree <- as.phylo(hclust(dist(t(mat_sig)), method="ward.D2"))
if("C2B" %in% phylo_tree$tip.label) phylo_tree <- root(phylo_tree, outgroup="C2B", resolve.root=TRUE)
grp_cols <- c("Culture_U2"="#1f77b4", "Primary_U2"="#ff7f0e", "Recurrent_U2"="#d62728")
tip_cols <- if("Classification" %in% colnames(meta)) grp_cols[as.character(meta[phylo_tree$tip.label, "Classification"])] else "black"
tip_cols[is.na(tip_cols)] <- "black"

tree_png <- paste0(output_prefix, "_Phylogenetic_Tree_mqc.png")
png(tree_png, width=10, height=8, units="in", res=300)
plot(phylo_tree, type="phylogram", tip.color=tip_cols, main="Phylogenetic Tree", cex=1.2, edge.width=2)
legend("topright", legend=names(grp_cols), fill=grp_cols, bty="n")
dev.off()
add_image_base64(tree_png, "Figure 1: Phylogenetic Tree (Ward's D2 Clustering)")

# --- PCA ---
pca <- prcomp(t(mat_sig))
var_pc <- round(summary(pca)$importance[2, 1:2]*100, 1)
pcaData <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Sample=rownames(meta), Class=meta$Classification)

centroids <- aggregate(cbind(PC1, PC2) ~ Class, data=pcaData, FUN=mean)
centroids <- centroids[match(levels(meta$Classification), centroids$Class), ]
centroids <- na.omit(centroids) 

arrow_data <- data.frame(
    x_start = centroids$PC1[-nrow(centroids)],
    y_start = centroids$PC2[-nrow(centroids)],
    x_end = centroids$PC1[-1],
    y_end = centroids$PC2[-1]
)

perm_res <- tryCatch({ adonis2(dist(t(mat_sig)) ~ Classification, data = meta, permutations = 999) }, error=function(e) NULL)
sub_txt <- if(!is.null(perm_res)) paste0("PERMANOVA: p=", perm_res$`Pr(>F)`[1], ", R2=", round(perm_res$R2[1]*100,1), "%") else ""

p_pca <- ggplot(pcaData, aes(x=PC1, y=PC2)) +
    geom_segment(data=arrow_data, aes(x=x_start, y=y_start, xend=x_end, yend=y_end),
                 arrow=arrow(length=unit(0.4,"cm"), type="closed"), color="grey50", linewidth=1.5, inherit.aes=FALSE) +
    geom_point(aes(fill=Class, shape=Class), size=7, color="black", stroke=0.8) +
    ggrepel::geom_text_repel(aes(label=Sample), size=4, box.padding=0.5, point.padding=0.3) +
    geom_label(data=centroids, aes(x=PC1, y=PC2, label=Class), fill="white", alpha=0.8, fontface="bold") +
    scale_fill_manual(values=grp_cols) +
    scale_shape_manual(values=c(21, 24, 22)) +
    labs(title="Evolutionary Trajectory (PCA)", subtitle=sub_txt, 
         x=paste0("PC1 (", var_pc[1], "%)"), y=paste0("PC2 (", var_pc[2], "%)")) +
    theme_bw(base_size=TEXT_BASE_SIZE) + theme(legend.position="bottom")

pca_png <- paste0(output_prefix, "_PCA_Trajectory_mqc.png")
ggsave(pca_png, p_pca, width=10, height=9, dpi=300)
add_image_base64(pca_png, "Figure 2: Evolutionary Trajectory. The arrows indicate the net movement of the transcriptome centroid from Culture -> Primary -> Recurrent.")

# ==============================================================================
# 3. SUBTYPE SCORING & RAW TABLE
# ==============================================================================
add_header("2. Subtype Scoring")
cat("LOG [3/8]: Scoring Subtypes...\n")
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
gsva_res <- suppressWarnings(gsva(mat_sym, sigs, method="gsva", kcdf="Gaussian", verbose=FALSE))
gsva_scaled <- t(scale(t(gsva_res)))

# ADD RAW DATA TABLE
raw_df <- as.data.frame(t(gsva_res))
raw_df <- tibble::rownames_to_column(raw_df, "Sample")
add_table(raw_df, "Raw Subtype Scores (GSVA)", scroll=TRUE)

# ==============================================================================
# 4. PLASTICITY (ENTROPY)
# ==============================================================================
add_header("3. Plasticity (Entropy Analysis)")
add_explainer("Shannon Entropy", "Measures the 'confusion' of a cell's identity. High entropy means the sample expresses a mix of all subtypes.")

cat("LOG [4/8]: Calculating Entropy...\n")
calc_entropy <- function(scores) {
    probs <- exp(scores) / sum(exp(scores))
    -sum(probs * log(probs))
}
meta$Plasticity <- apply(gsva_res, 2, calc_entropy)

plast_aov <- summary(aov(Plasticity ~ Classification, data=meta))
plast_p <- plast_aov[[1]][["Pr(>F)"]][1]

p_plast <- ggplot(meta, aes(x=Classification, y=Plasticity, fill=Classification)) +
    geom_boxplot(alpha=0.6) + geom_jitter(width=0.1, size=3) +
    scale_fill_manual(values=grp_cols) +
    labs(title="Subtype Plasticity", subtitle=paste0("ANOVA P = ", safe_format(plast_p))) +
    theme_bw(base_size=14)

plast_png <- paste0(output_prefix, "_Plasticity_Entropy_mqc.png")
ggsave(plast_png, p_plast, width=7, height=6, dpi=300)
add_image_base64(plast_png, "Figure 3: Plasticity Score. Higher entropy suggests cells are less differentiated.")

# ==============================================================================
# 5. HIGH-SENSITIVITY STATS (LIMMA & VARIANCE)
# ==============================================================================
add_header("4. Statistical Analysis (Limma & Variance)")
add_explainer("Limma (Empirical Bayes)", "Detects differential subtype activity.")
add_explainer("Bartlett's Test", "Tests if subtype DIVERSITY (variance) changes.")

cat("LOG [5/8]: Running Limma...\n")

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
trend_res <- limma_res %>% dplyr::filter(P.Value < 0.15) %>% 
    dplyr::select(Subtype, Comparison, logFC, P.Value, adj.P.Val) %>% arrange(P.Value)
trend_res$P.Value <- sapply(trend_res$P.Value, safe_format)
trend_res$adj.P.Val <- sapply(trend_res$adj.P.Val, safe_format)
add_table(trend_res, "Differential Trajectories (Limma P < 0.15)")

var_df <- data.frame()
for(sig in rownames(gsva_res)) {
    scores <- gsva_res[sig, ]
    bp <- tryCatch(bartlett.test(scores ~ meta$Classification)$p.value, error=function(e) NA)
    var_df <- rbind(var_df, data.frame(Subtype=sig, Bartlett_Variance_P=safe_format(bp)))
}
add_table(var_df, "Heterogeneity Tests (Variance Differences)")

# ==============================================================================
# 6. GLOBAL HEATMAP
# ==============================================================================
add_header("5. Subtype Landscape")
cat("LOG [6/8]: Generating Heatmap...\n")

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

ht_png <- paste0(output_prefix, "_Subtype_Heatmap_mqc.png")
png(ht_png, width=10, height=10, units="in", res=300)
draw(ht)
dev.off()
add_image_base64(ht_png, "Figure 4: Global Subtype Heatmap.")

# ==============================================================================
# 7. LLM SUMMARY & FINALIZE
# ==============================================================================
llm_text <- paste0(
    "SUMMARY FOR LLM:\n",
    "- Global Structure: PERMANOVA p-value = ", sub_txt, "\n",
    "- Plasticity (Entropy): ANOVA p-value = ", safe_format(plast_p), "\n",
    "- Top Differential Subtypes (Limma P<0.15): ", paste(trend_res$Subtype[1:3], collapse=", "), "\n",
    "- Heterogeneity Change (Bartlett P<0.05): ", paste(var_df$Subtype[as.numeric(var_df$Bartlett_Variance_P)<0.05], collapse=", "), "\n"
)
add_llm_summary(llm_text)

add_citations()
finish_html(output_prefix)
cat("\n=== PIPELINE COMPLETE ===\n")
