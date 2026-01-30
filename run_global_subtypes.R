#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# GLOBAL SUBTYPES ANALYSIS v13.0 (Robust Small-N Edition)
# ------------------------------------------------------------------------------

set.seed(12345)

suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(ape); library(ggrepel)
    library(EnsDb.Hsapiens.v86); library(ComplexHeatmap); library(circlize)
    library(GSVA); library(tidyr); library(tibble)
    library(limma); library(vegan); library(car)
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

safe_format <- function(x) {
    if (is.null(x) || length(x) == 0) return("NA")
    if (length(x) > 1) return(paste(x, collapse=", "))  # Handle vectors
    if (is.na(x)) return("NA")
    if (!is.numeric(x)) return(as.character(x))
    return(formatC(x, format="e", digits=2))
}

interpret_p <- function(p) {
    if (is.null(p) || length(p) == 0 || is.na(p)) return("NA")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    if (p < 0.10) return(".")
    return("ns")
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

cat("LOG [1/12]: Loading Data...\n")
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

# ==============================================================================
# 2. GLOBAL STRUCTURE
# ==============================================================================
cat("LOG [2/12]: Global Structure...\n")
top_var <- head(order(apply(mat_vst, 1, var), decreasing=TRUE), N_TOP_VARIABLE_GENES)
mat_sig <- mat_vst[top_var, ]

# PCA & PERMANOVA
perm_res <- tryCatch({
    adonis2(dist(t(mat_sig)) ~ Classification, data = meta, permutations = 999)
}, error=function(e) NULL)

perm_p <- if(!is.null(perm_res)) perm_res$`Pr(>F)`[1] else NA
perm_r2 <- if(!is.null(perm_res)) perm_res$R2[1] else NA

pca <- prcomp(t(mat_sig))
pcaData <- data.frame(
    PC1 = pca$x[,1], 
    PC2 = pca$x[,2], 
    Sample = rownames(meta), 
    Class = meta$Classification
)

n_pcs <- min(10, ncol(pca$x))
var_pc <- round(summary(pca)$importance[2, 1:n_pcs]*100, 1)

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

# Trajectory Plot
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

p_pca <- ggplot(pcaData, aes(x=PC1, y=PC2)) +
    {if(nrow(arrow_data) > 0) 
        geom_segment(data=arrow_data, 
                    aes(x=x_start, y=y_start, xend=x_end, yend=y_end), 
                    arrow=arrow(length=unit(0.4,"cm"), type="closed"), 
                    color="grey50", linewidth=1.2, inherit.aes=FALSE)
    } +
    geom_point(aes(fill=Class, shape=Class), size=6, color="black", stroke=0.8) +
    ggrepel::geom_text_repel(aes(label=Sample), size=4, box.padding=0.5) +
    scale_fill_manual(values=c("#1f77b4", "#ff7f0e", "#d62728")) + 
    scale_shape_manual(values=c(21, 24, 22)) +
    labs(title="Evolutionary Trajectory", 
         subtitle=paste0("PERMANOVA: R²=", round(perm_r2, 3), ", P=", safe_format(perm_p)), 
         x=paste0("PC1 (", var_pc[1], "%)"), 
         y=paste0("PC2 (", var_pc[2], "%)")) + 
    theme_bw(base_size=TEXT_BASE_SIZE)

ggsave(paste0(output_prefix, "_PCA_Trajectory_mqc.png"), p_pca, width=9, height=8)

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
cat("LOG [3/12]: Dual Scoring Methods (GSVA + Z-Score)...\n")

sigs <- list(
    "Verhaak_MES" = c("CHI3L1","CD44","VIM","RELB","STAT3"), 
    "Verhaak_CL" = c("EGFR","AKT2","NOTCH3","CCND2"),
    "Neftel_AC" = c("APOE","AQP4","CLU","S100B"), 
    "Neftel_MES" = c("CHI3L1","CD44","ANXA1","VIM"),
    "Garofano_MTC" = c("SLC45A1","GOT2","IDH2","SDHA","CS"), 
    "Garofano_GPM" = c("SLC2A1","HK2","PFKP","LDHA","GAPDH")
)

# Method 1: GSVA
gsva_res <- suppressWarnings(
    gsva(mat_sym, sigs, method="gsva", kcdf="Gaussian", verbose=FALSE)
)

# Method 2: Z-Score Averaging
mat_z <- t(scale(t(mat_sym)))
z_res <- matrix(NA, nrow=length(sigs), ncol=ncol(mat_z), 
                dimnames=list(names(sigs), colnames(mat_z)))

for(s in names(sigs)) {
    genes <- intersect(sigs[[s]], rownames(mat_z))
    if(length(genes) > 1) {
        z_res[s,] <- colMeans(mat_z[genes, , drop=FALSE])
    } else if(length(genes) == 1) {
        z_res[s,] <- mat_z[genes, ]
    }
}

# Method agreement
if(SCORING_METHOD == "both") {
    cor_methods <- cor(as.vector(gsva_res), as.vector(z_res), use="complete.obs")
    cat(sprintf("\n✓ GSVA vs Z-score correlation: r = %.3f\n", cor_methods))
    
    agreement_df <- data.frame(
        GSVA = as.vector(gsva_res),
        Zscore = as.vector(z_res),
        Signature = rep(rownames(gsva_res), each=ncol(gsva_res))
    )
    
    p_agreement <- ggplot(agreement_df, aes(x=GSVA, y=Zscore, color=Signature)) +
        geom_point(alpha=0.7, size=3) +
        geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
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
} else {
    final_scores <- gsva_res
    scoring_label <- "GSVA"
}

cat(sprintf("→ Using %s for downstream analysis\n", scoring_label))

# ==============================================================================
# 4. PLASTICITY (ENTROPY)
# ==============================================================================
cat("LOG [4/12]: Plasticity Analysis...\n")

calc_entropy <- function(s) { 
    p <- exp(s)/sum(exp(s))
    -sum(p*log(p), na.rm=TRUE)
}

meta$Plasticity <- apply(final_scores, 2, calc_entropy)

shapiro_p <- shapiro.test(meta$Plasticity)$p.value

if(shapiro_p < 0.05) {
    plast_test <- kruskal.test(Plasticity ~ Classification, data=meta)
    plast_p <- plast_test$p.value
    plast_method <- "Kruskal-Wallis"
} else {
    plast_aov <- aov(Plasticity ~ Classification, data=meta)
    plast_p <- summary(plast_aov)[[1]][["Pr(>F)"]][1]
    plast_method <- "ANOVA"
}

p_plast <- ggplot(meta, aes(x=Classification, y=Plasticity, fill=Classification)) +
    geom_boxplot(alpha=0.6, outlier.shape=NA) + 
    geom_jitter(width=0.1, size=3, alpha=0.7) +
    scale_fill_manual(values=c("#1f77b4", "#ff7f0e", "#d62728")) +
    labs(title="Cellular Plasticity (Shannon Entropy)", 
         subtitle=paste0(plast_method, " P=", safe_format(plast_p))) + 
    theme_bw(base_size=14) +
    theme(legend.position="none")

ggsave(paste0(output_prefix, "_Plasticity_Entropy_mqc.png"), p_plast, width=6, height=6)

# ==============================================================================
# 5. HEATMAPS & CORRELATION
# ==============================================================================
cat("LOG [5/12]: Heatmaps & Co-evolution...\n")

col_ann <- HeatmapAnnotation(
    Class = meta$Classification,
    col = list(Class = c(
        "Culture_U2" = "#1f77b4", 
        "Primary_U2" = "#ff7f0e", 
        "Recurrent_U2" = "#d62728"
    ))
)

ht <- Heatmap(
    final_scores, 
    name = "Score", 
    top_annotation = col_ann,
    col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), 
    cluster_columns = FALSE,
    column_order = order(meta$Classification)
)

png(paste0(output_prefix, "_Score_Heatmap_mqc.png"), 
    width=10, height=6, units="in", res=300)
draw(ht)
dev.off()

sig_cor <- cor(t(final_scores), method="pearson")

ht_cor <- Heatmap(
    sig_cor, 
    name = "Pearson\nCorr", 
    col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    column_title = "Signature Co-evolution"
)

png(paste0(output_prefix, "_Sig_Correlation_mqc.png"), 
    width=7, height=6, units="in", res=300)
draw(ht_cor)
dev.off()

# ==============================================================================
# 6. ROBUST LIMMA (WITH ARRAY WEIGHTS)
# ==============================================================================
cat("LOG [6/12]: Robust Stats with arrayWeights...\n")

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

# ==============================================================================
# 7. MULTI-FDR SENSITIVITY & VOLCANOES
# ==============================================================================
cat("LOG [7/12]: Multi-FDR Sensitivity Analysis...\n")

sig_counts <- data.frame(FDR_Threshold=FDR_THRESHOLDS, Label=FDR_LABELS)
for(c in all_coefs) sig_counts[[c]] <- 0

for(i in 1:length(FDR_THRESHOLDS)) {
    thresh <- FDR_THRESHOLDS[i]
    lbl <- FDR_LABELS[i]
    
    stats_out <- data.frame(Signature=rownames(final_scores))
    
    for(coef in all_coefs) {
        tt <- topTable(fit, coef=coef, number=Inf, adjust.method="BH")
        
        stats_out[[paste0(coef, "_logFC")]] <- round(tt[stats_out$Signature, "logFC"], 3)
        stats_out[[paste0(coef, "_adjP")]] <- sapply(stats_out$Signature, function(sig) {
            safe_format(tt[sig, "adj.P.Val"])
        })
        stats_out[[paste0(coef, "_Sig")]] <- ifelse(tt[stats_out$Signature, "adj.P.Val"] < thresh, "Yes", "No")
        
        sig_counts[i, coef] <- sum(tt[, "adj.P.Val"] < thresh, na.rm=TRUE)
        
        volc_data <- data.frame(
            Signature = rownames(tt), 
            logFC = tt$logFC, 
            negLogP = -log10(tt$adj.P.Val), 
            Significant = tt$adj.P.Val < thresh
        )
        
        p_vol <- ggplot(volc_data, aes(x=logFC, y=negLogP, color=Significant)) +
            geom_point(size=4, alpha=0.7) + 
            geom_text_repel(data=subset(volc_data, Significant), 
                           aes(label=Signature), size=3.5) +
            geom_hline(yintercept=-log10(thresh), linetype="dashed", color="red") +
            geom_vline(xintercept=0, linetype="dashed", color="grey50") +
            scale_color_manual(values=c("grey60", "red")) + 
            theme_bw(base_size=12) +
            labs(title=paste("Volcano:", coef), subtitle=lbl, 
                 x="log2 Fold Change", y="-log10(Adj P)")
        
        ggsave(paste0(output_prefix, "_Volcano_", gsub("_", "", coef), "_FDR", 
                     gsub("[^0-9]", "", as.character(thresh)), "_mqc.png"), 
               p_vol, width=7, height=6)
    }
    
    tf <- paste0(dirname(output_prefix), "/fdr_", 
                gsub("[^0-9]", "", as.character(thresh)), "_stats_mqc.tsv")
    cat(paste0("# id: 'stats_", gsub("[^0-9]", "", as.character(thresh)), 
               "'\n# section_name: 'Stats: ", lbl, 
               "'\n# plot_type: 'table'\n"), file=tf)
    write.table(stats_out, file=tf, sep="\t", quote=FALSE, row.names=FALSE, append=TRUE)
}

# ==============================================================================
# 8. TRAJECTORY PLOTS
# ==============================================================================
cat("LOG [8/12]: Signature Trajectories...\n")

for(sig in rownames(final_scores)) {
    df <- data.frame(
        Score = final_scores[sig, ], 
        Class = meta$Classification
    )
    
    summ <- df %>% 
        group_by(Class) %>% 
        summarise(Mean = mean(Score)) %>% 
        mutate(Num = as.numeric(Class))
    
    p <- ggplot(df, aes(x=as.numeric(Class), y=Score)) +
        geom_line(data=summ, aes(x=Num, y=Mean, group=1), 
                 color="darkred", linewidth=1.2, alpha=0.7) +
        geom_point(aes(fill=Class), size=4, shape=21, alpha=0.8) + 
        scale_fill_manual(values=c("#1f77b4", "#ff7f0e", "#d62728")) +
        scale_x_continuous(breaks=1:3, labels=levels(meta$Classification)) +
        theme_bw(base_size=12) + 
        labs(title=sig, y="Score", x="Stage") +
        theme(legend.position="none")
    
    ggsave(paste0(output_prefix, "_Trajectory_", sig, "_mqc.png"), p, width=6, height=4)
}

# ==============================================================================
# 9. SENSITIVITY SUMMARY
# ==============================================================================
cat("LOG [9/12]: FDR Sensitivity Summary...\n")

sig_counts_long <- sig_counts %>%
    pivot_longer(cols=all_of(all_coefs), names_to="Contrast", values_to="N_Significant") %>%
    mutate(Contrast = factor(Contrast, levels=all_coefs))

p_sensitivity <- ggplot(sig_counts_long, aes(x=Label, y=N_Significant, fill=Contrast)) +
    geom_bar(stat="identity", position="dodge", alpha=0.8) +
    geom_text(aes(label=N_Significant), position=position_dodge(0.9), vjust=-0.5, size=3) +
    scale_fill_manual(values=c("#1f77b4", "#ff7f0e", "#d62728")) +
    labs(title="FDR Sensitivity",
         subtitle=paste0(nrow(final_scores), " signatures tested"),
         x="FDR Threshold", y="N Significant") +
    theme_bw(base_size=12) +
    theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "bottom") +
    ylim(0, nrow(final_scores) + 0.5)

ggsave(paste0(output_prefix, "_FDR_Sensitivity_mqc.png"), p_sensitivity, width=10, height=7)

summary_file <- paste0(dirname(output_prefix), "/fdr_sensitivity_summary_mqc.tsv")
cat("# id: 'fdr_sensitivity'\n# section_name: 'FDR Sensitivity'\n# plot_type: 'table'\n", 
    file=summary_file)
write.table(sig_counts, file=summary_file, sep="\t", quote=FALSE, row.names=FALSE, append=TRUE)

# ==============================================================================
# 10. EXPORT DATA
# ==============================================================================
cat("LOG [10/12]: Exporting Data...\n")

write.csv(t(final_scores), paste0(output_prefix, "_Scores.csv"))
write.csv(meta, paste0(output_prefix, "_Metadata.csv"))
write.csv(weight_df, paste0(output_prefix, "_SampleWeights.csv"), row.names=FALSE)

# ==============================================================================
# 11. HTML SUMMARY
# ==============================================================================
cat("LOG [11/12]: Generating HTML Summary...\n")

summary_html <- paste0(dirname(output_prefix), "/analysis_summary_mqc.html")
sink(summary_html)

# Create simple text summaries
group_summary <- paste(names(group_sizes), as.vector(group_sizes), sep=": ", collapse="\n")

cat("
<div style='background-color:#f8f9fa; padding:15px; border:1px solid #ddd; margin-bottom:20px;'>
    <h3>🔬 Analysis Summary (Robust Small-N Edition)</h3>
    
    <h4>Metadata</h4>
    <pre>", group_summary, "</pre>
    
    <h4>Scoring Method</h4>
    <p><strong>", scoring_label, "</strong></p>
    ", if(SCORING_METHOD == "both") paste0("<p>Method correlation: r = ", round(cor_methods, 3), "</p>") else "", "
    
    <h4>Global Structure</h4>
    <ul>
        <li>PERMANOVA: R²=", round(perm_r2, 3), ", P=", safe_format(perm_p), " ", interpret_p(perm_p), "</li>
        <li>Plasticity (", plast_method, "): P=", safe_format(plast_p), " ", interpret_p(plast_p), "</li>
    </ul>
    
    <h4>Sample Quality Weights</h4>
    <table style='border-collapse: collapse; margin:10px 0;'>
        <tr style='background:#ddd;'><th style='padding:5px; border:1px solid #999;'>Sample</th>
            <th style='padding:5px; border:1px solid #999;'>Weight</th>
            <th style='padding:5px; border:1px solid #999;'>Group</th></tr>
")

for(j in order(weight_df$Weight)) {
    cat(sprintf("        <tr><td style='padding:5px; border:1px solid #ccc;'>%s</td>
            <td style='padding:5px; border:1px solid #ccc;'>%.3f</td>
            <td style='padding:5px; border:1px solid #ccc;'>%s</td></tr>\n",
            weight_df$Sample[j], weight_df$Weight[j], weight_df$Group[j]))
}

cat("    </table>
    
    <h4>PCA Drivers</h4>
    <pre style='font-size:10px;'>", pc_drivers, "</pre>
    
    <h4>FDR Sensitivity Results</h4>
    <table style='border-collapse: collapse; width:100%;'>
        <tr style='background:#ddd;'><th style='padding:5px; border:1px solid #999;'>Threshold</th>")

for(coef in all_coefs) {
    cat(sprintf("<th style='padding:5px; border:1px solid #999;'>%s</th>", gsub("_", " ", coef)))
}
cat("</tr>\n")

for(i in 1:nrow(sig_counts)) {
    cat(sprintf("        <tr><td style='padding:5px; border:1px solid #ccc;'><strong>%s</strong></td>", 
                sig_counts$Label[i]))
    for(coef in all_coefs) {
        cat(sprintf("<td style='padding:5px; border:1px solid #ccc; text-align:center;'>%d</td>", 
                   sig_counts[i, coef]))
    }
    cat("</tr>\n")
}

cat("    </table>
</div>

<div style='font-size:11px; color:#666; margin-top:20px;'>
    <h4>Methods</h4>
    <ul>
        <li><strong>Scoring:</strong> ", scoring_label, "</li>
        <li><strong>Statistical Model:</strong> limma with arrayWeights</li>
        <li><strong>FDR Correction:</strong> Benjamini-Hochberg</li>
        <li><strong>Software:</strong> R ", as.character(R.version.string), 
        ", limma ", as.character(packageVersion("limma")), "</li>
    </ul>
</div>
")

sink()

# ==============================================================================
# 12. SESSION INFO
# ==============================================================================
cat("LOG [12/12]: Saving Session Info...\n")
writeLines(capture.output(sessionInfo()), paste0(dirname(output_prefix), "/sessionInfo.txt"))

cat("\n=== PIPELINE COMPLETE ===\n")
cat("✓ Dual scoring with agreement metrics\n")
cat("✓ Sample quality weights calculated\n")
cat("✓ Signature co-evolution matrix\n")
cat("\nFDR Summary:\n")
print(sig_counts)
