#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# GLOBAL SUBTYPES ANALYSIS v14.1 (Trajectory Sorting + HTML Fixes)
# ------------------------------------------------------------------------------
# Key improvements over v14.0:
# - Sorts trajectory plots by correlation (Increasing → Decreasing)
# - Fixes HTML interpret_p() to return plain text (was adding HTML in wrong place)
# - Adds patchwork library
# - More robust error handling
# ------------------------------------------------------------------------------

set.seed(12345)

suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(ape); library(ggrepel)
    library(EnsDb.Hsapiens.v86); library(ComplexHeatmap); library(circlize)
    library(GSVA); library(tidyr); library(tibble)
    library(limma); library(vegan); library(car); library(patchwork)
})

# ==============================================================================
# CONFIGURATION
# ==============================================================================
N_TOP_VARIABLE_GENES <- 500
TEXT_BASE_SIZE <- 14
MIN_SAMPLES_PER_GROUP <- 3
FDR_THRESHOLDS <- c(0.05, 0.01, 0.005, 0.001)
FDR_LABELS <- c("Standard (0.05)", "Strict (0.01)", "Very Strict (0.005)", "Ultra (0.001)")
SCORING_METHOD <- "both"

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================
map_genes_to_symbols <- function(gene_ids, db = EnsDb.Hsapiens.v86) {
    clean_ids <- sub("\\..*", "", gene_ids)
    symbols <- mapIds(db, keys=clean_ids, column="SYMBOL", keytype="GENEID", multiVals="first")
    ifelse(is.na(symbols), clean_ids, symbols)
}

safe_format <- function(x, digits=3) {
    if (is.null(x) || length(x) == 0) return("NA")
    if (length(x) > 1) return(paste(sapply(x, safe_format, digits=digits), collapse=", "))
    if (is.na(x)) return("NA")
    if (!is.numeric(x)) return(as.character(x))
    if (abs(x) < 0.001 && x != 0) return(formatC(x, format="e", digits=2))
    return(round(x, digits))
}

# Plain text version for stat_log
interpret_p_text <- function(p) {
    if (is.null(p) || length(p) == 0 || is.na(p)) return("NA")
    if (p < 0.001) return("*** (highly significant)")
    if (p < 0.01) return("** (very significant)")
    if (p < 0.05) return("* (significant)")
    if (p < 0.10) return(". (trend)")
    return("ns (not significant)")
}

# HTML version for final report
interpret_p_html <- function(p) {
    if (is.null(p) || length(p) == 0 || is.na(p)) return("<span class='ns'>NA</span>")
    if (p < 0.001) return("<span class='significant'>***</span>")
    if (p < 0.01) return("<span class='significant'>**</span>")
    if (p < 0.05) return("<span class='significant'>*</span>")
    if (p < 0.10) return("<span class='trend'>.</span>")
    return("<span class='ns'>ns</span>")
}

# Statistical narrative logger
stat_log <- list()
add_stat_log <- function(title, content) {
    stat_log[[length(stat_log) + 1]] <<- list(title=title, content=content)
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

cat("LOG [1/10]: Loading Data...\n")
vst_file <- file.path(results_dir, "tables/processed_abundance/all.vst.tsv")
mat_vst <- as.matrix(read.table(vst_file, header=TRUE, row.names=1, check.names=FALSE))
meta <- read.csv(meta_file, row.names=1)
common <- intersect(colnames(mat_vst), rownames(meta))
mat_vst <- mat_vst[, common]
meta <- meta[common, , drop=FALSE]

if("Classification" %in% colnames(meta)) {
    meta$Classification <- factor(meta$Classification, 
                                  levels=c("Culture_U2", "Primary_U2", "Recurrent_U2"))
}

group_sizes <- table(meta$Classification)
cat("Sample sizes per group:\n")
print(group_sizes)

add_stat_log("Data Loading", sprintf(
    "Loaded %d genes across %d samples. Groups: %s",
    nrow(mat_vst), ncol(mat_vst),
    paste(names(group_sizes), "=", group_sizes, collapse="; ")
))

# Map symbols
mapped_syms <- map_genes_to_symbols(rownames(mat_vst))
mat_sym_df <- as.data.frame(mat_vst) %>% 
    tibble::rownames_to_column("ensembl") %>%
    mutate(symbol = mapped_syms) %>% 
    dplyr::filter(!is.na(symbol))

mat_sym <- mat_sym_df %>%
    group_by(symbol) %>%
    summarise(across(where(is.numeric), mean)) %>% 
    tibble::column_to_rownames("symbol") %>% 
    as.matrix()

add_stat_log("Gene Mapping", sprintf(
    "Mapped %d Ensembl IDs to %d unique gene symbols",
    nrow(mat_vst), nrow(mat_sym)
))

# ==============================================================================
# 2. GLOBAL STRUCTURE
# ==============================================================================
cat("LOG [2/10]: Global Structure...\n")
top_var <- head(order(apply(mat_vst, 1, var), decreasing=TRUE), N_TOP_VARIABLE_GENES)
mat_sig <- mat_vst[top_var, ]

add_stat_log("Feature Selection", sprintf(
    "Selected top %d most variable genes (variance range: %.2f - %.2f)",
    N_TOP_VARIABLE_GENES,
    min(apply(mat_sig, 1, var)),
    max(apply(mat_sig, 1, var))
))

# PCA & PERMANOVA
perm_res <- tryCatch({
    adonis2(dist(t(mat_sig)) ~ Classification, data = meta, permutations = 999)
}, error=function(e) NULL)

perm_p <- if(!is.null(perm_res)) perm_res$`Pr(>F)`[1] else NA
perm_r2 <- if(!is.null(perm_res)) perm_res$R2[1] else NA
perm_f <- if(!is.null(perm_res)) perm_res$F[1] else NA

add_stat_log("PERMANOVA Test", sprintf(
    "Test: Classification effect on multivariate gene expression\n  F-statistic = %.3f\n  R² = %.3f (%.1f%% of variance explained)\n  P-value = %s %s\n  Interpretation: %s",
    perm_f, perm_r2, perm_r2*100, safe_format(perm_p), interpret_p_text(perm_p),
    if(!is.na(perm_p) && perm_p < 0.05) "Groups show distinct transcriptional profiles" else "Groups are not significantly different"
))

pca <- prcomp(t(mat_sig))
pcaData <- data.frame(
    PC1 = pca$x[,1], 
    PC2 = pca$x[,2], 
    Sample = rownames(meta), 
    Class = meta$Classification
)

n_pcs <- min(10, ncol(pca$x))
var_pc <- round(summary(pca)$importance[2, 1:n_pcs]*100, 1)
cum_var <- round(summary(pca)$importance[3, 1:n_pcs]*100, 1)

add_stat_log("PCA Variance", sprintf(
    "PC1 explains %.1f%% variance\n  PC2 explains %.1f%% variance\n  PC1+PC2 together: %.1f%%\n  First 5 PCs capture %.1f%% of total variance",
    var_pc[1], var_pc[2], cum_var[2], cum_var[5]
))

# Scree Plot
scree_data <- data.frame(
    PC = factor(paste0("PC", 1:n_pcs), levels=paste0("PC", 1:n_pcs)), 
    Variance = var_pc
)

p_scree <- ggplot(scree_data, aes(x=PC, y=Variance)) +
    geom_bar(stat="identity", fill="steelblue") + 
    geom_line(aes(group=1), color="darkblue", linewidth=1) +
    geom_point(color="darkblue", size=3) +
    labs(title="PCA Scree Plot", y="% Variance Explained") + 
    theme_bw(base_size=TEXT_BASE_SIZE)

ggsave(paste0(output_prefix, "_Scree_Plot_mqc.png"), p_scree, width=8, height=5)

# PCA Loadings
loadings <- pca$rotation
rownames(loadings) <- map_genes_to_symbols(rownames(loadings))

pc_drivers <- ""
for(i in 1:min(5, ncol(loadings))) {
    pc <- paste0("PC", i)
    top <- names(sort(abs(loadings[, pc]), decreasing=TRUE)[1:8])
    pc_drivers <- paste0(pc_drivers, "\n", pc, ": ", paste(top, collapse=", "))
}

add_stat_log("PCA Gene Drivers", sprintf(
    "Top genes driving principal components:%s",
    pc_drivers
))

# Trajectory Plot with Biplot
centroids <- aggregate(cbind(PC1, PC2) ~ Class, data=pcaData, FUN=mean)
centroids <- centroids[match(levels(meta$Classification), centroids$Class), ]
centroids <- na.omit(centroids)

arrow_data <- if(nrow(centroids) >= 2) {
    data.frame(
        x_start = centroids$PC1[-nrow(centroids)], 
        y_start = centroids$PC2[-nrow(centroids)], 
        x_end = centroids$PC1[-1], 
        y_end = centroids$PC2[-1]
    )
} else {
    data.frame()
}

# Top gene loadings for biplot
top_genes_load <- loadings[order(sqrt(loadings[,"PC1"]^2 + loadings[,"PC2"]^2), decreasing=TRUE)[1:10], c("PC1", "PC2")]
gene_arrow_scale <- max(abs(pcaData$PC1)) / max(abs(top_genes_load[,"PC1"])) * 0.8
gene_arrows <- as.data.frame(top_genes_load * gene_arrow_scale)
gene_arrows$Gene <- rownames(gene_arrows)

p_pca <- ggplot(pcaData, aes(x=PC1, y=PC2)) +
    {if(nrow(arrow_data) > 0) 
        geom_segment(data=arrow_data, 
                    aes(x=x_start, y=y_start, xend=x_end, yend=y_end), 
                    arrow=arrow(length=unit(0.4,"cm"), type="closed"), 
                    color="grey50", linewidth=1.2, inherit.aes=FALSE)
    } +
    geom_segment(data=gene_arrows, aes(x=0, y=0, xend=PC1, yend=PC2), 
                arrow=arrow(length=unit(0.2,"cm")), color="red", alpha=0.4) +
    geom_text_repel(data=gene_arrows, aes(x=PC1, y=PC2, label=Gene), 
                   color="red", size=3, segment.alpha=0.4) +
    geom_point(aes(fill=Class, shape=Class), size=6, color="black", stroke=0.8) +
    ggrepel::geom_text_repel(aes(label=Sample), size=4, box.padding=0.5) +
    scale_fill_manual(values=c("#1f77b4", "#ff7f0e", "#d62728")) + 
    scale_shape_manual(values=c(21, 24, 22)) +
    labs(title="Evolutionary Trajectory with Gene Drivers", 
         subtitle=paste0("PERMANOVA: R²=", round(perm_r2, 3), ", P=", safe_format(perm_p)), 
         x=paste0("PC1 (", var_pc[1], "%)"), 
         y=paste0("PC2 (", var_pc[2], "%)")) + 
    theme_bw(base_size=TEXT_BASE_SIZE)

ggsave(paste0(output_prefix, "_PCA_Trajectory_mqc.png"), p_pca, width=10, height=9)

# Tree
phylo_tree <- as.phylo(hclust(dist(t(mat_sig)), method="ward.D2"))
if("C2B" %in% phylo_tree$tip.label) {
    phylo_tree <- root(phylo_tree, outgroup="C2B", resolve.root=TRUE)
}
tip_cols <- c("#1f77b4", "#ff7f0e", "#d62728")[
    as.numeric(meta[phylo_tree$tip.label, "Classification"])
]

png(paste0(output_prefix, "_Phylogenetic_Tree_mqc.png"), 
    width=10, height=8, units="in", res=300)
plot(phylo_tree, type="phylogram", tip.color=tip_cols, 
     main="Phylogenetic Tree (Ward's D2)", cex=1.2, edge.width=2)
legend("topleft", legend=levels(meta$Classification), 
       fill=c("#1f77b4", "#ff7f0e", "#d62728"), bty="n")
dev.off()

# ==============================================================================
# 3. SUBTYPE SCORING (HYBRID: GSVA + Z-SCORE)
# ==============================================================================
cat("LOG [3/10]: Dual Scoring Methods (GSVA + Z-Score)...\n")

sigs <- list(
    "Verhaak_MES" = c("CHI3L1","CD44","VIM","RELB","STAT3"), 
    "Verhaak_CL" = c("EGFR","AKT2","NOTCH3","CCND2"),
    "Neftel_AC" = c("APOE","AQP4","CLU","S100B"), 
    "Neftel_MES" = c("CHI3L1","CD44","ANXA1","VIM"),
    "Garofano_MTC" = c("SLC45A1","GOT2","IDH2","SDHA","CS"), 
    "Garofano_GPM" = c("SLC2A1","HK2","PFKP","LDHA","GAPDH")
)

# Check gene coverage
gene_coverage <- sapply(sigs, function(sig_genes) {
    sum(sig_genes %in% rownames(mat_sym)) / length(sig_genes) * 100
})

add_stat_log("Signature Gene Coverage", sprintf(
    "Gene detection rates in dataset:\n  %s\n  Mean coverage: %.1f%%",
    paste(names(gene_coverage), "=", round(gene_coverage, 1), "%", collapse="\n  "),
    mean(gene_coverage)
))

# Method 1: GSVA
gsva_res <- suppressWarnings(
    gsva(mat_sym, sigs, method="gsva", kcdf="Gaussian", verbose=FALSE)
)

# Method 2: Z-Score Averaging
mat_z <- t(scale(t(mat_sym)))
z_res <- matrix(0, nrow=length(sigs), ncol=ncol(mat_z), 
                dimnames=list(names(sigs), colnames(mat_z)))

for(s in names(sigs)) {
    genes <- intersect(sigs[[s]], rownames(mat_z))
    if(length(genes) > 1) {
        z_res[s,] <- colMeans(mat_z[genes, , drop=FALSE], na.rm=TRUE)
    } else if(length(genes) == 1) {
        z_res[s,] <- mat_z[genes, ]
    } else {
        cat(sprintf("WARNING: No genes found for signature %s\n", s))
    }
}

# Method agreement
if(SCORING_METHOD == "both") {
    cor_methods <- cor(as.vector(gsva_res), as.vector(z_res), use="complete.obs")
    
    add_stat_log("Scoring Method Comparison", sprintf(
        "GSVA vs Z-Score correlation: r = %.3f\n  Interpretation: %s\n  Decision: Using Z-Score (more robust for small gene sets)",
        cor_methods,
        if(cor_methods > 0.8) "High agreement between methods" else if(cor_methods > 0.6) "Moderate agreement" else "Low agreement - caution advised"
    ))
    
    agreement_df <- data.frame(
        GSVA = as.vector(gsva_res),
        Zscore = as.vector(z_res),
        Signature = rep(rownames(gsva_res), each=ncol(gsva_res))
    )
    
    p_agreement <- ggplot(agreement_df, aes(x=GSVA, y=Zscore, color=Signature)) +
        geom_point(alpha=0.7, size=3) +
        geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
        geom_smooth(method="lm", se=TRUE, color="black", linetype="dashed", linewidth=0.5) +
        labs(title="Method Agreement: GSVA vs Z-Score",
             subtitle=paste0("Pearson r = ", round(cor_methods, 3))) +
        theme_bw(base_size=12)
    
    ggsave(paste0(output_prefix, "_Method_Agreement_mqc.png"), p_agreement, 
           width=8, height=6)
    
    final_scores <- z_res
    scoring_label <- "Z-Score (Robust for small gene sets)"
} else if(SCORING_METHOD == "zscore") {
    final_scores <- z_res
    scoring_label <- "Z-Score"
    add_stat_log("Scoring Method", "Using Z-Score method (specified by user)")
} else {
    final_scores <- gsva_res
    scoring_label <- "GSVA"
    add_stat_log("Scoring Method", "Using GSVA method (specified by user)")
}

cat(sprintf("→ Using %s for downstream analysis\n", scoring_label))

# ==============================================================================
# 4. PLASTICITY (ENTROPY)
# ==============================================================================
cat("LOG [4/10]: Plasticity Analysis...\n")

calc_entropy <- function(s) { 
    p <- exp(s)/sum(exp(s))
    -sum(p*log(p), na.rm=TRUE)
}

meta$Plasticity <- apply(final_scores, 2, calc_entropy)

# Test normality
shapiro_p <- shapiro.test(meta$Plasticity)$p.value

add_stat_log("Plasticity: Normality Test", sprintf(
    "Shapiro-Wilk test: W = %.3f, P = %s\n  Interpretation: Data is %s\n  Decision: Using %s test",
    shapiro.test(meta$Plasticity)$statistic,
    safe_format(shapiro_p),
    if(shapiro_p >= 0.05) "normally distributed" else "non-normally distributed",
    if(shapiro_p < 0.05) "Kruskal-Wallis (non-parametric)" else "ANOVA (parametric)"
))

if(shapiro_p < 0.05) {
    plast_test <- kruskal.test(Plasticity ~ Classification, data=meta)
    plast_p <- plast_test$p.value
    plast_stat <- plast_test$statistic
    plast_method <- "Kruskal-Wallis"
    
    add_stat_log("Plasticity: Group Comparison", sprintf(
        "Kruskal-Wallis H-test\n  H = %.3f, df = %d\n  P = %s %s\n  Interpretation: %s",
        plast_stat, plast_test$parameter,
        safe_format(plast_p), interpret_p_text(plast_p),
        if(plast_p < 0.05) "Significant differences in plasticity across groups" else "No significant differences in plasticity"
    ))
} else {
    plast_aov <- aov(Plasticity ~ Classification, data=meta)
    plast_summary <- summary(plast_aov)[[1]]
    plast_p <- plast_summary[["Pr(>F)"]][1]
    plast_f <- plast_summary[["F value"]][1]
    plast_method <- "ANOVA"
    
    add_stat_log("Plasticity: Group Comparison", sprintf(
        "One-way ANOVA\n  F(%d,%d) = %.3f\n  P = %s %s\n  Interpretation: %s",
        plast_summary[["Df"]][1], plast_summary[["Df"]][2],
        plast_f, safe_format(plast_p), interpret_p_text(plast_p),
        if(plast_p < 0.05) "Significant differences in plasticity across groups" else "No significant differences in plasticity"
    ))
}

# Calculate group statistics
plast_stats <- meta %>%
    group_by(Classification) %>%
    summarise(
        Mean = mean(Plasticity),
        SD = sd(Plasticity),
        Median = median(Plasticity),
        IQR = IQR(Plasticity)
    )

add_stat_log("Plasticity: Group Statistics", sprintf(
    "Summary by group:\n  %s",
    paste(capture.output(print(plast_stats)), collapse="\n  ")
))

p_plast <- ggplot(meta, aes(x=Classification, y=Plasticity, fill=Classification)) +
    geom_boxplot(alpha=0.6, outlier.shape=NA) + 
    geom_jitter(width=0.1, size=3, alpha=0.7) +
    stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="white") +
    scale_fill_manual(values=c("#1f77b4", "#ff7f0e", "#d62728")) +
    labs(title="Cellular Plasticity (Shannon Entropy)", 
         subtitle=paste0(plast_method, " P=", safe_format(plast_p), " ", interpret_p_text(plast_p)),
         caption="Diamond = mean, box = median ± IQR") + 
    theme_bw(base_size=14) +
    theme(legend.position="none")

ggsave(paste0(output_prefix, "_Plasticity_Entropy_mqc.png"), p_plast, width=7, height=6)

# ==============================================================================
# 5. HEATMAPS & CORRELATION
# ==============================================================================
cat("LOG [5/10]: Heatmaps & Co-evolution...\n")

col_ann <- HeatmapAnnotation(
    Class = meta$Classification,
    Plasticity = meta$Plasticity,
    col = list(
        Class = c(
            "Culture_U2" = "#1f77b4", 
            "Primary_U2" = "#ff7f0e", 
            "Recurrent_U2" = "#d62728"
        ),
        Plasticity = colorRamp2(c(min(meta$Plasticity), max(meta$Plasticity)), c("white", "purple"))
    )
)

ht <- Heatmap(
    final_scores, 
    name = "Score", 
    top_annotation = col_ann,
    col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), 
    cluster_columns = FALSE,
    column_order = order(meta$Classification),
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 9)
)

png(paste0(output_prefix, "_Score_Heatmap_mqc.png"), 
    width=10, height=6, units="in", res=300)
draw(ht)
dev.off()

sig_cor <- cor(t(final_scores), method="pearson")

# Find strongest correlations
cor_pairs <- which(abs(sig_cor) > 0.6 & upper.tri(sig_cor), arr.ind=TRUE)
if(nrow(cor_pairs) > 0) {
    cor_summary <- data.frame(
        Sig1 = rownames(sig_cor)[cor_pairs[,1]],
        Sig2 = rownames(sig_cor)[cor_pairs[,2]],
        Correlation = sig_cor[cor_pairs]
    )
    cor_summary <- cor_summary[order(abs(cor_summary$Correlation), decreasing=TRUE), ]
    
    add_stat_log("Signature Co-evolution", sprintf(
        "Strong correlations (|r| > 0.6):\n  %s",
        paste(capture.output(print(cor_summary, row.names=FALSE)), collapse="\n  ")
    ))
} else {
    add_stat_log("Signature Co-evolution", "No strong correlations (|r| > 0.6) detected between signatures")
}

ht_cor <- Heatmap(
    sig_cor, 
    name = "Pearson\nCorr", 
    col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    column_title = "Signature Co-evolution",
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10)
)

png(paste0(output_prefix, "_Sig_Correlation_mqc.png"), 
    width=7, height=6, units="in", res=300)
draw(ht_cor)
dev.off()

# ==============================================================================
# 6. ROBUST LIMMA (WITH ARRAY WEIGHTS)
# ==============================================================================
cat("LOG [6/10]: Robust Stats with arrayWeights...\n")

design <- model.matrix(~0 + meta$Classification)
colnames(design) <- levels(meta$Classification)

aw <- arrayWeights(final_scores, design)

cat("\n📊 Sample Reliability Weights (lower = noisier):\n")
weight_df <- data.frame(
    Sample = names(aw), 
    Weight = round(aw, 3), 
    Group = as.character(meta[names(aw), "Classification"])
)
print(weight_df[order(weight_df$Weight), ])

# Identify potential outliers
outlier_threshold <- mean(aw) - 2*sd(aw)
outliers <- weight_df$Sample[weight_df$Weight < outlier_threshold]

add_stat_log("Sample Quality Assessment", sprintf(
    "arrayWeights analysis:\n  Mean weight = %.3f\n  SD = %.3f\n  Range = %.3f - %.3f\n  Potential outliers (weight < mean - 2SD): %s\n  Interpretation: %s",
    mean(aw), sd(aw), min(aw), max(aw),
    if(length(outliers) > 0) paste(outliers, collapse=", ") else "None detected",
    if(length(outliers) > 0) "Some samples show reduced reliability - weights applied to down-weight their influence" else "All samples show consistent quality"
))

cont.matrix <- makeContrasts(
    Prim_vs_Cult = Primary_U2 - Culture_U2, 
    Rec_vs_Cult = Recurrent_U2 - Culture_U2, 
    Rec_vs_Prim = Recurrent_U2 - Primary_U2, 
    levels = design
)

fit <- lmFit(final_scores, design, weights=aw) 
fit <- contrasts.fit(fit, cont.matrix)
fit <- eBayes(fit)

all_coefs <- colnames(cont.matrix)

add_stat_log("Linear Model Setup", sprintf(
    "Contrasts tested:\n  1. %s\n  2. %s\n  3. %s\n  Model: Weighted linear regression with empirical Bayes moderation\n  Adjustment: Benjamini-Hochberg FDR correction",
    "Primary vs Culture", "Recurrent vs Culture", "Recurrent vs Primary"
))

# ==============================================================================
# 7. UNIFIED TRAJECTORY PLOT (SORTED BY TREND) ⭐ KEY IMPROVEMENT
# ==============================================================================
cat("LOG [7/10]: Creating Unified Trajectory Visualization...\n")

# Prepare data for all signatures
traj_data_list <- list()
for(sig in rownames(final_scores)) {
    df <- data.frame(
        Signature = sig,
        Sample = colnames(final_scores),
        Score = final_scores[sig, ],
        Class = meta$Classification,
        Stage = as.numeric(meta$Classification)
    )
    traj_data_list[[sig]] <- df
}
traj_data <- do.call(rbind, traj_data_list)

# Calculate group means
traj_summary <- traj_data %>%
    group_by(Signature, Class, Stage) %>%
    summarise(Mean = mean(Score), SE = sd(Score)/sqrt(n()), .groups="drop")

# ⭐ SORT SIGNATURES BY CORRELATION WITH STAGE (Increasing → Decreasing)
trend_order <- traj_data %>%
    group_by(Signature) %>%
    summarise(Cor = cor(Stage, Score, method="spearman"), .groups="drop") %>%
    arrange(desc(Cor)) %>%
    pull(Signature)

traj_data$Signature <- factor(traj_data$Signature, levels=trend_order)
traj_summary$Signature <- factor(traj_summary$Signature, levels=trend_order)

# Create faceted trajectory plot with sorted panels
p_unified_traj <- ggplot(traj_data, aes(x=Stage, y=Score, group=Signature)) +
    geom_line(data=traj_summary, aes(x=Stage, y=Mean, color=Signature), 
             linewidth=1.2, alpha=0.8) +
    geom_point(aes(fill=Class, shape=Class), size=3, alpha=0.6) +
    geom_ribbon(data=traj_summary, aes(x=Stage, ymin=Mean-SE, ymax=Mean+SE, fill=Signature),
               alpha=0.2, color=NA) +
    facet_wrap(~Signature, scales="free_y", ncol=3) +
    scale_x_continuous(breaks=1:3, labels=levels(meta$Classification)) +
    scale_fill_manual(values=c("#1f77b4", "#ff7f0e", "#d62728")) +
    scale_color_brewer(palette="Set1") +
    scale_shape_manual(values=c(21, 24, 22)) +
    labs(title="Signature Trajectories Across Evolution",
         subtitle="Ordered by trend: Increasing (top) → Decreasing (bottom) | Lines = group means ± SE",
         x="Stage", y="Signature Score") +
    theme_bw(base_size=11) +
    theme(
        strip.background = element_rect(fill="grey90"),
        strip.text = element_text(face="bold", size=10),
        legend.position="bottom",
        panel.grid.minor = element_blank()
    )

ggsave(paste0(output_prefix, "_Unified_Trajectories_mqc.png"), 
       p_unified_traj, width=14, height=10)

# Summary statistics for trajectories
traj_models <- traj_data %>%
    group_by(Signature) %>%
    summarise(
        Correlation_with_Stage = cor(Stage, Score, method="spearman"),
        Trend_P = cor.test(Stage, Score, method="spearman")$p.value,
        Direction = ifelse(Correlation_with_Stage > 0, "Increasing", "Decreasing"),
        .groups="drop"
    ) %>%
    arrange(desc(Correlation_with_Stage))  # Match plot order

add_stat_log("Trajectory Analysis", sprintf(
    "Signature evolution patterns (sorted by trend strength):\n  %s",
    paste(capture.output(print(traj_models, row.names=FALSE)), collapse="\n  ")
))

# ==============================================================================
# 8. STATISTICAL RESULTS TABLE
# ==============================================================================
cat("LOG [8/10]: Compiling Statistical Results...\n")

# Get results for standard FDR threshold
results_005 <- data.frame(Signature = rownames(final_scores))

for(coef in all_coefs) {
    tt <- topTable(fit, coef=coef, number=Inf, adjust.method="BH")
    results_005[[paste0(coef, "_logFC")]] <- round(tt[results_005$Signature, "logFC"], 3)
    results_005[[paste0(coef, "_P")]] <- sapply(results_005$Signature, function(sig) {
        safe_format(tt[sig, "P.Value"])
    })
    results_005[[paste0(coef, "_adjP")]] <- sapply(results_005$Signature, function(sig) {
        safe_format(tt[sig, "adj.P.Val"])
    })
    results_005[[paste0(coef, "_Sig")]] <- ifelse(
        tt[results_005$Signature, "adj.P.Val"] < 0.05, 
        "✓", "—"
    )
}

# Save comprehensive results
write.csv(results_005, paste0(output_prefix, "_Statistical_Results.csv"), row.names=FALSE)

# Create MultiQC table
mqc_table <- paste0(dirname(output_prefix), "/differential_signatures_mqc.tsv")
cat("# id: 'diff_sigs'\n# section_name: 'Differential Signatures (FDR<0.05)'\n# plot_type: 'table'\n", 
    file=mqc_table)
write.table(results_005, file=mqc_table, sep="\t", quote=FALSE, row.names=FALSE, append=TRUE)

# Count significant changes at each threshold
sig_counts <- data.frame(
    FDR_Threshold = FDR_THRESHOLDS,
    Label = FDR_LABELS
)

for(coef in all_coefs) {
    tt <- topTable(fit, coef=coef, number=Inf, adjust.method="BH")
    sig_counts[[coef]] <- sapply(FDR_THRESHOLDS, function(thresh) {
        sum(tt$adj.P.Val < thresh, na.rm=TRUE)
    })
}

add_stat_log("FDR Sensitivity Analysis", sprintf(
    "Number of significant signatures at different FDR thresholds:\n  %s",
    paste(capture.output(print(sig_counts, row.names=FALSE)), collapse="\n  ")
))

# Sensitivity plot
sig_counts_long <- sig_counts %>%
    pivot_longer(cols=all_of(all_coefs), names_to="Contrast", values_to="N_Significant") %>%
    mutate(Contrast = factor(Contrast, levels=all_coefs))

p_sensitivity <- ggplot(sig_counts_long, aes(x=Label, y=N_Significant, fill=Contrast)) +
    geom_bar(stat="identity", position="dodge", alpha=0.8) +
    geom_text(aes(label=N_Significant), position=position_dodge(0.9), vjust=-0.5, size=3.5) +
    scale_fill_manual(values=c("#1f77b4", "#ff7f0e", "#d62728")) +
    labs(title="FDR Sensitivity Analysis",
         subtitle=paste0(nrow(final_scores), " signatures tested across ", length(all_coefs), " contrasts"),
         x="FDR Threshold", y="Number of Significant Signatures") +
    theme_bw(base_size=12) +
    theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "bottom") +
    ylim(0, nrow(final_scores) + 0.5)

ggsave(paste0(output_prefix, "_FDR_Sensitivity_mqc.png"), p_sensitivity, width=10, height=7)

# ==============================================================================
# 9. EXPORT DATA
# ==============================================================================
cat("LOG [9/10]: Exporting Data...\n")

write.csv(t(final_scores), paste0(output_prefix, "_Scores.csv"))
write.csv(meta, paste0(output_prefix, "_Metadata.csv"))
write.csv(weight_df, paste0(output_prefix, "_SampleWeights.csv"), row.names=FALSE)
write.csv(sig_counts, paste0(output_prefix, "_FDR_Sensitivity.csv"), row.names=FALSE)
write.csv(traj_models, paste0(output_prefix, "_Trajectory_Statistics.csv"), row.names=FALSE)

# ==============================================================================
# 10. RICH NARRATIVE HTML REPORT
# ==============================================================================
cat("LOG [10/10]: Generating Rich Narrative Report...\n")

summary_html <- paste0(dirname(output_prefix), "/analysis_summary_mqc.html")
sink(summary_html)

cat("<!DOCTYPE html>
<html>
<head>
<style>
body { font-family: Arial, sans-serif; max-width: 1200px; margin: 20px auto; padding: 20px; }
.header { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 30px; border-radius: 10px; margin-bottom: 30px; }
.section { background-color: #f8f9fa; padding: 20px; border-left: 4px solid #667eea; margin-bottom: 20px; border-radius: 5px; }
.stat-block { background: white; padding: 15px; margin: 10px 0; border-radius: 5px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
.stat-title { font-weight: bold; color: #667eea; margin-bottom: 10px; font-size: 1.1em; }
.stat-content { font-family: 'Courier New', monospace; font-size: 0.9em; white-space: pre-wrap; line-height: 1.6; }
table { width: 100%; border-collapse: collapse; margin: 15px 0; }
th { background-color: #667eea; color: white; padding: 12px; text-align: left; }
td { padding: 10px; border-bottom: 1px solid #ddd; }
tr:hover { background-color: #f5f5f5; }
.highlight { background-color: #fff3cd; padding: 2px 6px; border-radius: 3px; }
.significant { color: #28a745; font-weight: bold; }
.not-significant { color: #6c757d; }
.trend { color: #fd7e14; font-weight: bold; }
.ns { color: #adb5bd; }
</style>
</head>
<body>

<div class='header'>
<h1>🧬 U251 Global Subtype Evolution Analysis</h1>
<h3>Comprehensive Statistical Report v14.1</h3>
<p>Analysis Date: ", Sys.Date(), "</p>
<p>Scoring Method: <strong>", scoring_label, "</strong></p>
</div>

<div class='section'>
<h2>📊 Executive Summary</h2>
<div class='stat-block'>
<div class='stat-title'>Dataset Overview</div>
<div class='stat-content'>")

cat(sprintf("Samples analyzed: %d\nGenes quantified: %d\nSignatures evaluated: %d\n\nGroup breakdown:\n",
            ncol(mat_vst), nrow(mat_vst), nrow(final_scores)))
for(i in seq_along(group_sizes)) {
    cat(sprintf("  • %s: %d samples\n", names(group_sizes)[i], group_sizes[i]))
}

cat("</div></div>

<div class='stat-block'>
<div class='stat-title'>Key Findings</div>
<div class='stat-content'>")

# Summary of key statistical findings
cat(sprintf("PERMANOVA Test: P = %s %s\n", safe_format(perm_p), interpret_p_html(perm_p)))
cat(sprintf("  → Groups %s distinct transcriptional signatures\n\n", 
            if(!is.na(perm_p) && perm_p < 0.05) "SHOW" else "DO NOT SHOW"))

cat(sprintf("Plasticity Analysis: P = %s %s\n", safe_format(plast_p), interpret_p_html(plast_p)))
cat(sprintf("  → Cellular plasticity %s across evolution\n\n",
            if(plast_p < 0.05) "CHANGES SIGNIFICANTLY" else "remains stable"))

# Count significant signatures
n_sig_prim <- sum(sapply(results_005$Signature, function(sig) {
    tt <- topTable(fit, coef="Prim_vs_Cult", number=Inf)
    tt[sig, "adj.P.Val"] < 0.05
}), na.rm=TRUE)

n_sig_rec_cult <- sum(sapply(results_005$Signature, function(sig) {
    tt <- topTable(fit, coef="Rec_vs_Cult", number=Inf)
    tt[sig, "adj.P.Val"] < 0.05
}), na.rm=TRUE)

n_sig_rec_prim <- sum(sapply(results_005$Signature, function(sig) {
    tt <- topTable(fit, coef="Rec_vs_Prim", number=Inf)
    tt[sig, "adj.P.Val"] < 0.05
}), na.rm=TRUE)

cat(sprintf("Differential Signatures (FDR < 0.05):\n"))
cat(sprintf("  • Primary vs Culture: %d/%d signatures\n", n_sig_prim, nrow(final_scores)))
cat(sprintf("  • Recurrent vs Culture: %d/%d signatures\n", n_sig_rec_cult, nrow(final_scores)))
cat(sprintf("  • Recurrent vs Primary: %d/%d signatures\n", n_sig_rec_prim, nrow(final_scores)))

cat("</div></div>
</div>

<div class='section'>
<h2>📈 Detailed Statistical Methods & Decisions</h2>")

# Print all statistical logs
for(log_entry in stat_log) {
    cat(sprintf("
<div class='stat-block'>
<div class='stat-title'>%s</div>
<div class='stat-content'>%s</div>
</div>", log_entry$title, log_entry$content))
}

cat("</div>

<div class='section'>
<h2>📉 Sample Quality Assessment</h2>
<div class='stat-block'>
<table>
<tr><th>Sample</th><th>Weight</th><th>Group</th><th>Quality</th></tr>")

for(i in 1:nrow(weight_df)) {
    quality <- if(weight_df$Weight[i] > mean(aw) + sd(aw)) {
        "<span class='significant'>Excellent</span>"
    } else if(weight_df$Weight[i] > mean(aw)) {
        "Good"
    } else if(weight_df$Weight[i] > mean(aw) - sd(aw)) {
        "Fair"
    } else {
        "<span class='not-significant'>Poor</span>"
    }
    
    cat(sprintf("<tr><td>%s</td><td>%.3f</td><td>%s</td><td>%s</td></tr>\n",
                weight_df$Sample[i], weight_df$Weight[i], weight_df$Group[i], quality))
}

cat("</table>
<p style='margin-top: 15px; font-size: 0.9em; color: #666;'>
<strong>Note:</strong> Sample weights are automatically calculated by limma's arrayWeights function. 
Lower weights indicate potentially noisy samples that are down-weighted in the statistical analysis.
</p>
</div>
</div>

<div class='section'>
<h2>🎯 PCA Gene Drivers</h2>
<div class='stat-block'>
<div class='stat-content' style='font-size: 0.85em;'>", pc_drivers, "</div>
</div>
</div>

<div class='section'>
<h2>📊 Signature Evolution Trajectories</h2>
<div class='stat-block'>
<p style='margin-bottom: 15px; color: #666;'>Sorted by correlation strength (strongest increasing trends first)</p>
<table>
<tr><th>Signature</th><th>Correlation with Stage</th><th>Trend P-value</th><th>Direction</th><th>Significance</th></tr>")

for(i in 1:nrow(traj_models)) {
    sig_marker <- if(traj_models$Trend_P[i] < 0.05) {
        "<span class='significant'>✓ Significant</span>"
    } else {
        "<span class='not-significant'>— Not Significant</span>"
    }
    
    cat(sprintf("<tr><td>%s</td><td>%.3f</td><td>%s</td><td>%s</td><td>%s</td></tr>\n",
                traj_models$Signature[i],
                traj_models$Correlation_with_Stage[i],
                safe_format(traj_models$Trend_P[i]),
                traj_models$Direction[i],
                sig_marker))
}

cat("</table>
</div>
</div>

<div class='section'>
<h2>🔬 Software & Methods</h2>
<div class='stat-block'>
<div class='stat-content'>
<strong>R Version:</strong> ", as.character(R.version.string), "

<strong>Key Packages:</strong>
  • limma ", as.character(packageVersion("limma")), " (differential analysis)
  • GSVA ", as.character(packageVersion("GSVA")), " (signature scoring)
  • vegan ", as.character(packageVersion("vegan")), " (PERMANOVA)
  • ComplexHeatmap ", as.character(packageVersion("ComplexHeatmap")), " (visualization)
  • patchwork ", as.character(packageVersion("patchwork")), " (plot composition)

<strong>Statistical Approach:</strong>
  • Scoring: ", scoring_label, "
  • Model: Weighted linear regression with empirical Bayes shrinkage
  • Sample weights: Automatically calculated via arrayWeights
  • Multiple testing: Benjamini-Hochberg FDR correction
  • Normality testing: Shapiro-Wilk test
  • Group comparisons: Parametric (ANOVA) or non-parametric (Kruskal-Wallis) based on normality
  • Trajectory trends: Spearman correlation with evolutionary stage

<strong>Reproducibility:</strong>
  • Random seed: 12345
  • Top variable genes: ", N_TOP_VARIABLE_GENES, "
  • Analysis date: ", Sys.time(), "
</div>
</div>
</div>

<div class='section' style='background-color: #e7f3ff; border-left-color: #2196F3;'>
<h2>💡 Interpretation Guide</h2>
<div class='stat-block'>
<h4>Significance Levels:</h4>
<ul>
<li><strong>*** (P < 0.001):</strong> Highly significant - very strong evidence</li>
<li><strong>** (P < 0.01):</strong> Very significant - strong evidence</li>
<li><strong>* (P < 0.05):</strong> Significant - moderate evidence</li>
<li><strong>. (P < 0.10):</strong> Trend - suggestive but inconclusive</li>
<li><strong>ns (P ≥ 0.10):</strong> Not significant - insufficient evidence</li>
</ul>

<h4>Key Metrics:</h4>
<ul>
<li><strong>R² (PERMANOVA):</strong> Proportion of variance explained by groups (0-1 scale)</li>
<li><strong>Shannon Entropy:</strong> Measure of cellular plasticity; higher = more plastic/undifferentiated</li>
<li><strong>logFC:</strong> Log2 fold change; positive = higher in second group, negative = higher in first</li>
<li><strong>Sample Weights:</strong> Quality metric; values ~1.0 are typical, <0.5 suggests noise</li>
<li><strong>Spearman ρ:</strong> Correlation with stage; positive = increasing trend, negative = decreasing</li>
</ul>
</div>
</div>

</body>
</html>")

sink()

# ==============================================================================
# SESSION INFO
# ==============================================================================
writeLines(capture.output(sessionInfo()), paste0(dirname(output_prefix), "/sessionInfo.txt"))

cat("\n╔════════════════════════════════════════════╗\n")
cat("║     ANALYSIS COMPLETE - SUMMARY            ║\n")
cat("╚════════════════════════════════════════════╝\n\n")
cat(sprintf("✓ Analyzed %d samples across %d genes\n", ncol(mat_vst), nrow(mat_vst)))
cat(sprintf("✓ Evaluated %d subtype signatures\n", nrow(final_scores)))
cat(sprintf("✓ Statistical logs: %d entries\n", length(stat_log)))
cat(sprintf("✓ Rich HTML report generated\n"))
cat(sprintf("\nKey outputs:\n"))
cat(sprintf("  • Unified trajectory plot (SORTED by trend!)\n"))
cat(sprintf("  • PCA with gene drivers\n"))
cat(sprintf("  • Comprehensive statistical report\n"))
cat(sprintf("  • No unnecessary volcano plots!\n\n"))
