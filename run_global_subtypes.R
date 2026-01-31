#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# GLOBAL SUBTYPES ANALYSIS v16.2 (DUAL-WEIGHTING + TRAJECTORY + OMNIBUS EDITION)
# ------------------------------------------------------------------------------
# NEW IN v16.2:
# 1. Fixed trajectory plot rendering (na.rm=TRUE in summarise) 
# 2. Dual scoring: GSVA + Z-Score with consensus reporting
# 3. Enhanced LLM summary with PCA, plasticity, correlations
# 4. All v16.1 features preserved (dual weighting, JT tests, etc.)
# ------------------------------------------------------------------------------

set.seed(12345)

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(dplyr)
    library(ape)
    library(ggrepel)
    library(EnsDb.Hsapiens.v86)
    library(ComplexHeatmap)
    library(circlize)
    library(GSVA)
    library(tidyr)
    library(tibble)
    library(limma)
    library(vegan)
    library(car)
    library(patchwork)
    library(RColorBrewer)
    library(clinfun)  # For Jonckheere-Terpstra test
})

# ==============================================================================
# CLI ARGUMENTS
# ==============================================================================
option_list <- list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="Path to VST abundance table", metavar="file"),
    make_option(c("-o", "--out"), type="character", default="results/subtypes/global",
                help="Output prefix", metavar="str"),
    make_option(c("-m", "--meta"), type="character", default=NULL,
                help="Path to metadata CSV", metavar="file"),
    make_option(c("-s", "--signatures"), type="character", default=NULL,
                help="Optional: Path to custom signatures (GMT or CSV)", metavar="file"),
    make_option(c("--n_top_genes"), type="integer", default=500,
                help="Number of top variable genes for PCA [default=%default]", metavar="int"),
    make_option(c("--scoring"), type="character", default="both",
                help="Scoring method: 'zscore', 'gsva', or 'both' [default=%default]", metavar="str")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$meta)) {
    print_help(opt_parser)
    stop("Input and Metadata files are required.", call.=FALSE)
}

# ==============================================================================
# CONFIGURATION
# ==============================================================================
N_TOP_VAR <- opt$n_top_genes
SCORING_METHOD <- opt$scoring
FDR_THRESHOLDS <- c(0.05, 0.01, 0.005, 0.001)
FDR_LABELS <- c("Standard (0.05)", "Strict (0.01)", "Very Strict (0.005)", "Ultra (0.001)")

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================
safe_format <- function(x, digits=3) {
    if (is.null(x) || length(x) == 0 || is.na(x)) return("NA")
    if (length(x) > 1) return(paste(sapply(x, safe_format, digits=digits), collapse=", "))
    if (!is.numeric(x)) return(as.character(x))
    if (abs(x) < 0.001 && x != 0) return(formatC(x, format="e", digits=2))
    return(round(x, digits))
}

interpret_p <- function(p) {
    if (is.null(p) || length(p) == 0 || is.na(p)) return("NA")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    if (p < 0.10) return(".")
    return("ns")
}

theme_publication <- function(base_size=14) {
    theme_bw(base_size=base_size) +
        theme(
            panel.grid.minor = element_blank(),
            plot.title = element_text(face="bold", size=rel(1.2)),
            plot.subtitle = element_text(color="grey40", size=rel(0.9)),
            legend.position = "bottom",
            strip.background = element_rect(fill="grey95", color=NA),
            strip.text = element_text(face="bold")
        )
}

stat_log <- list()
add_stat_log <- function(title, content) {
    stat_log[[length(stat_log) + 1]] <<- list(title=title, content=content)
}

# ==============================================================================
# 1. DATA LOADING
# ==============================================================================
cat("LOG [1/11]: Loading Data...\n")

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)

mat_vst <- as.matrix(read.table(opt$input, header=TRUE, row.names=1, check.names=FALSE))
meta <- read.csv(opt$meta, row.names=1)

common <- intersect(colnames(mat_vst), rownames(meta))
if(length(common) < 3) stop("Error: <3 matching samples between VST and Metadata.")
mat_vst <- mat_vst[, common]
meta <- meta[common, , drop=FALSE]

if(!"Classification" %in% colnames(meta)) {
    stop("Metadata must contain a 'Classification' column.")
}

default_groups <- c("Culture_U2", "Primary_U2", "Recurrent_U2")
if(all(default_groups %in% unique(meta$Classification))) {
    meta$Classification <- factor(meta$Classification, levels=default_groups)
} else {
    meta$Classification <- as.factor(meta$Classification)
}

group_sizes <- table(meta$Classification)
cat("Sample sizes per group:\n")
print(group_sizes)

add_stat_log("Data Loading", sprintf(
    "Loaded %d genes across %d samples. Groups: %s",
    nrow(mat_vst), ncol(mat_vst),
    paste(names(group_sizes), "=", group_sizes, collapse="; ")
))

default_colors <- c("Culture_U2" = "#1f77b4", "Primary_U2" = "#ff7f0e", "Recurrent_U2" = "#d62728")
default_shapes <- c("Culture_U2" = 21, "Primary_U2" = 24, "Recurrent_U2" = 22)
current_groups <- levels(meta$Classification)

if(all(current_groups %in% names(default_colors))) {
    GROUP_COLORS <- default_colors[current_groups]
    GROUP_SHAPES <- default_shapes[current_groups]
    cat("  ✓ Using default U251 color scheme\n")
} else {
    cat("  ! Custom groups detected. Generating dynamic palette...\n")
    n_groups <- length(current_groups)
    pal <- brewer.pal(max(3, n_groups), "Set1")[1:n_groups]
    GROUP_COLORS <- setNames(pal, current_groups)
    GROUP_SHAPES <- setNames(rep(c(21, 24, 22, 23, 25), length.out=n_groups), current_groups)
}

sample_id <- rownames(mat_vst)[1]
if (grepl("^ENSG", sample_id)) {
    cat("  → Detected Ensembl IDs. Mapping to gene symbols...\n")
    clean_ids <- sub("\\..*", "", rownames(mat_vst))
    symbols <- mapIds(EnsDb.Hsapiens.v86, keys=clean_ids, column="SYMBOL",
                     keytype="GENEID", multiVals="first")

    mat_sym_df <- as.data.frame(mat_vst) %>%
        tibble::rownames_to_column("ensembl") %>%
        mutate(symbol = ifelse(is.na(symbols), clean_ids, symbols)) %>%
        dplyr::filter(!is.na(symbol)) %>%
        group_by(symbol) %>%
        summarise(across(where(is.numeric), mean)) %>%
        tibble::column_to_rownames("symbol")
    mat_sym <- as.matrix(mat_sym_df)

    add_stat_log("Gene Mapping", sprintf(
        "Mapped %d Ensembl IDs to %d unique gene symbols",
        nrow(mat_vst), nrow(mat_sym)
    ))
} else {
    cat("  → Gene symbols detected. Skipping Ensembl mapping.\n")
    mat_sym <- mat_vst
    add_stat_log("Gene Mapping", "Input already in gene symbol format")
}

# ==============================================================================
# 2. GLOBAL STRUCTURE (PCA + PERMANOVA + BIPLOT)
# ==============================================================================
cat("LOG [2/11]: Global Structure Analysis...\n")

top_var <- head(order(apply(mat_sym, 1, var), decreasing=TRUE), N_TOP_VAR)
mat_sig <- mat_sym[top_var, ]

add_stat_log("Feature Selection", sprintf(
    "Selected top %d most variable genes (variance range: %.2f - %.2f)",
    N_TOP_VAR, min(apply(mat_sig, 1, var)), max(apply(mat_sig, 1, var))
))

perm_res <- tryCatch({
    adonis2(dist(t(mat_sig)) ~ Classification, data = meta, permutations = 999)
}, error=function(e) NULL)

perm_p <- if(!is.null(perm_res)) perm_res$`Pr(>F)`[1] else NA
perm_r2 <- if(!is.null(perm_res)) perm_res$R2[1] else NA
perm_f <- if(!is.null(perm_res)) perm_res$F[1] else NA

add_stat_log("PERMANOVA Test", sprintf(
    "F = %.3f, R² = %.3f (%.1f%% variance explained by groups), P = %s %s\n  Interpretation: Groups %s significantly different transcriptional profiles",
    perm_f, perm_r2, perm_r2*100, safe_format(perm_p), interpret_p(perm_p),
    if(!is.na(perm_p) && perm_p < 0.05) "have" else "do not have"
))

pca <- prcomp(t(mat_sig))
n_samples <- ncol(mat_sig)
max_pcs <- min(10, n_samples - 1, ncol(pca$x))

pca_summary <- summary(pca)
n_pcs_to_show <- min(max_pcs, ncol(pca_summary$importance))

var_pc <- round(pca_summary$importance[2, 1:n_pcs_to_show]*100, 1)
cum_var <- round(pca_summary$importance[3, 1:n_pcs_to_show]*100, 1)
pcaData <- data.frame(pca$x[,1:2], Sample=rownames(meta), Class=meta$Classification)

add_stat_log("PCA Variance", sprintf(
    "PC1 = %.1f%%, PC2 = %.1f%%, PC1+PC2 = %.1f%%\n  First %d PCs explain %.1f%% total variance",
    var_pc[1], var_pc[2], cum_var[2], min(5, n_pcs_to_show), cum_var[min(5, n_pcs_to_show)]
))

loadings <- pca$rotation
pc_drivers_text <- ""
pc_drivers_list <- list()

for(i in 1:min(5, ncol(loadings))) {
    pc <- paste0("PC", i)
    top <- names(sort(abs(loadings[, i]), decreasing=TRUE)[1:8])
    pc_drivers_list[[pc]] <- top
    pc_drivers_text <- paste0(pc_drivers_text, "\n", pc, ": ", paste(top, collapse=", "))
}

add_stat_log("PCA Gene Drivers", sprintf("Top genes driving each PC:%s", pc_drivers_text))

scree_df <- data.frame(
    PC = factor(paste0("PC", 1:n_pcs_to_show), levels=paste0("PC", 1:n_pcs_to_show)),
    Var = var_pc
)
p_scree <- ggplot(scree_df, aes(x=PC, y=Var)) +
    geom_bar(stat="identity", fill="steelblue", alpha=0.8) +
    geom_line(aes(group=1), color="darkblue", linewidth=1) +
    geom_point(color="darkblue", size=2) +
    labs(title="Variance Explained", y="% Variance", x=NULL) +
    theme_publication(base_size=12)

centroids <- aggregate(cbind(PC1, PC2) ~ Class, data=pcaData, FUN=mean)
centroids <- centroids[match(levels(meta$Classification), centroids$Class), ]
centroids <- na.omit(centroids)

arrow_data <- if(nrow(centroids) >= 2) {
    data.frame(
        x = centroids$PC1[-nrow(centroids)],
        y = centroids$PC2[-nrow(centroids)],
        xend = centroids$PC1[-1],
        yend = centroids$PC2[-1]
    )
} else {
    data.frame()
}

top_genes_load <- loadings[order(sqrt(loadings[,"PC1"]^2 + loadings[,"PC2"]^2), decreasing=TRUE)[1:10], c("PC1", "PC2")]
gene_arrow_scale <- max(abs(pcaData$PC1)) / max(abs(top_genes_load[,"PC1"])) * 0.8
gene_arrows <- as.data.frame(top_genes_load * gene_arrow_scale)
gene_arrows$Gene <- rownames(gene_arrows)

p_pca <- ggplot(pcaData, aes(x=PC1, y=PC2)) +
    {if(nrow(arrow_data) > 0)
        geom_segment(data=arrow_data, aes(x=x, y=y, xend=xend, yend=yend),
                    arrow=arrow(length=unit(0.4,"cm"), type="closed"),
                    color="grey50", linewidth=1.2, inherit.aes=FALSE)
    } +
    geom_segment(data=gene_arrows, aes(x=0, y=0, xend=PC1, yend=PC2),
                arrow=arrow(length=unit(0.2,"cm")), color="red", alpha=0.4, inherit.aes=FALSE) +
    geom_text_repel(data=gene_arrows, aes(x=PC1, y=PC2, label=Gene),
                   color="red", size=3, segment.alpha=0.4) +
    geom_point(aes(fill=Class, shape=Class), size=6, color="black", stroke=0.8) +
    geom_text_repel(aes(label=Sample), size=3.5, box.padding=0.5) +
    scale_fill_manual(values=GROUP_COLORS) +
    scale_shape_manual(values=GROUP_SHAPES) +
    labs(title="Evolutionary Trajectory with Gene Drivers",
         subtitle=sprintf("PERMANOVA: R²=%.3f, P=%s | Red arrows = gene loadings", perm_r2, safe_format(perm_p)),
         x=paste0("PC1 (", var_pc[1], "%)"),
         y=paste0("PC2 (", var_pc[2], "%)")) +
    theme_publication()

p_composite <- (p_pca | p_scree) + plot_layout(widths = c(2.5, 1))
ggsave(paste0(opt$out, "_Global_Structure_mqc.png"), p_composite, width=13, height=6)

phylo_tree <- as.phylo(hclust(dist(t(mat_sig)), method="ward.D2"))
tip_cols <- GROUP_COLORS[as.character(meta[phylo_tree$tip.label, "Classification"])]

png(paste0(opt$out, "_Phylogenetic_Tree_mqc.png"), width=10, height=8, units="in", res=300)
plot(phylo_tree, type="phylogram", tip.color=tip_cols,
     main="Phylogenetic Tree (Ward's D2)", cex=1.2, edge.width=2)
legend("topleft", legend=levels(meta$Classification), fill=GROUP_COLORS, bty="n")
dev.off()

# ==============================================================================
# 3. SUBTYPE SCORING (DUAL METHOD: GSVA + Z-SCORE)
# ==============================================================================
cat("LOG [3/11]: Scoring Signatures (DUAL: GSVA + Z-Score)...\n")

sigs <- list()

if (!is.null(opt$signatures)) {
    cat("  → Parsing external signature file:", opt$signatures, "\n")
    
    if (grepl("\\.gmt$", opt$signatures, ignore.case = TRUE)) {
        lines <- readLines(opt$signatures)
        for (line in lines) {
            parts <- strsplit(line, "\t")[[1]]
            if (length(parts) > 2) {
                sig_name <- parts[1]
                genes <- parts[3:length(parts)]
                sigs[[sig_name]] <- genes[genes != ""]
            }
        }
        cat(sprintf("  ✓ Loaded %d signatures from GMT file\n", length(sigs)))
    } 
    else if (grepl("\\.(csv|tsv)$", opt$signatures, ignore.case = TRUE)) {
        sep_char <- if(grepl("\\.csv$", opt$signatures)) "," else "\t"
        sig_df <- read.table(opt$signatures, header=TRUE, sep=sep_char, stringsAsFactors=FALSE)
        if(ncol(sig_df) >= 2) {
            sigs <- split(sig_df[[1]], sig_df[[2]])
            cat(sprintf("  ✓ Loaded %d signatures from CSV/TSV file\n", length(sigs)))
        } else {
            stop("CSV/TSV must have at least 2 columns (Gene, Signature)")
        }
    } else {
        stop("Unknown signature file format. Use .gmt, .csv, or .tsv")
    }
    
    add_stat_log("Signature Source", sprintf("Loaded %d custom signatures from: %s", length(sigs), basename(opt$signatures)))
} else {
    cat("  → Using COMPLETE built-in GBM signatures (Verhaak + Neftel + Garofano)\n")
    sigs <- list(
        "Verhaak_Classical" = c("EGFR", "NES", "NOTCH3", "JAG1", "HES5", "AKT2", "FGFR2", 
                                "CCND2", "CDK4", "RB1", "SOX2", "SFRP2"),
        "Verhaak_Mesenchymal" = c("CHI3L1", "CD44", "VIM", "RELB", "STAT3", "TRADD", "CASP1",
                                  "TNFRSF1A", "IKBKB", "NFKB1", "MET", "FOSL2", "TIMP1",
                                  "S100A4", "FN1", "COL1A2", "COL3A1", "SERPINE1"),
        "Verhaak_Proneural" = c("OLIG2", "DLL3", "ASCL1", "TCF12", "DCX", "SOX11", "PDGFRA",
                                "NKX2-2", "ERBB3", "NCAM1", "NRCAM", "NKX6-1", "PAX3", "NEUROG1"),
        "Verhaak_Neural" = c("NEFL", "GABRA1", "SYT1", "SLC12A5", "GRIA2", "SYN1", "GABBR1",
                             "SNAP25", "ATP1A3", "STX1A", "GABRG2", "NRXN1", "NLGN1"),
        "Neftel_AC" = c("APOE", "AQP4", "CLU", "S100B", "SLC1A2", "SLC1A3", "GFAP",
                        "ALDOC", "GLUL", "GJA1", "ATP1A2", "ATP1B2", "FGFR3", "CD9",
                        "MT3", "AGT", "PTN", "CST3", "SPARC", "SPARCL1"),
        "Neftel_OPC" = c("PDGFRA", "OLIG1", "OLIG2", "CSPG4", "SOX10", "PCDH15", "BCAN",
                         "VCAN", "TNR", "LHFPL3", "GPR17", "MEGF11", "APOD", "COL9A1",
                         "LUZP2", "PLP1", "MAG", "MBP", "MOG", "MOBP"),
        "Neftel_NPC" = c("DCX", "DLL3", "ASCL1", "NEUROG2", "STMN2", "SOX11", "SOX4",
                         "TUBB3", "ELAVL2", "ELAVL3", "ELAVL4", "MAP2", "MAPT", "NCAM1",
                         "NRCAM", "ROBO2", "L1CAM", "GAP43", "NEFM", "PRPH"),
        "Neftel_MES" = c("CHI3L1", "CD44", "ANXA1", "VIM", "S100A4", "S100A6", "S100A10",
                         "LGALS1", "LGALS3", "SERPINE1", "TIMP1", "FN1", "COL1A1", "COL1A2",
                         "COL3A1", "COL5A1", "LAMC1", "ITGA5", "ITGAV", "MMP2"),
        "Garofano_MTC" = c("CS", "ACO2", "IDH2", "IDH3A", "OGDH", "SUCLA2", "SDHA", "SDHB",
                           "FH", "MDH2", "NDUFS1", "NDUFS2", "NDUFV1", "UQCRC1", "UQCRC2",
                           "COX4I1", "COX5A", "COX7A2", "ATP5F1A", "ATP5F1B"),
        "Garofano_GPM" = c("SLC2A1", "SLC2A3", "HK1", "HK2", "GPI", "PFKP", "PFKL", "PFKM",
                           "ALDOA", "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "PKM",
                           "LDHA", "LDHB", "SLC16A1", "SLC16A3", "G6PD", "PGD"),
        "Garofano_NEU" = c("GAD1", "GAD2", "SLC1A1", "SLC1A2", "SLC1A3", "GLUL", "GLS",
                           "SYP", "SNAP25", "SYT1", "VAMP2", "STX1A", "ATP1A1", "ATP1A2",
                           "ATP1A3", "SLC12A5", "KCNJ10", "KCNJ11")
    )
    
    add_stat_log("Signature Source", sprintf(
        "Using complete built-in signatures: %d total (Verhaak: 4, Neftel: 4, Garofano: 3)",
        length(sigs)
    ))
}

gene_coverage <- sapply(sigs, function(sig_genes) {
    sum(sig_genes %in% rownames(mat_sym)) / length(sig_genes) * 100
})

add_stat_log("Signature Gene Coverage", sprintf(
    "Mean coverage: %.1f%%\n  %s",
    mean(gene_coverage),
    paste(names(gene_coverage), "=", round(gene_coverage, 1), "%", collapse="\n  ")
))

# Calculate GSVA
cat("  → Computing GSVA scores...\n")
gsva_res <- suppressWarnings(gsva(mat_sym, sigs, method="gsva", kcdf="Gaussian", verbose=FALSE))

# Calculate Z-Scores
cat("  → Computing Z-Scores...\n")
mat_z <- t(scale(t(mat_sym)))
z_res <- matrix(0, nrow=length(sigs), ncol=ncol(mat_z), dimnames=list(names(sigs), colnames(mat_z)))

for(s in names(sigs)) {
    genes <- intersect(sigs[[s]], rownames(mat_z))
    if(length(genes) > 1) {
        z_res[s,] <- colMeans(mat_z[genes, , drop=FALSE], na.rm=TRUE)
    } else if(length(genes) == 1) {
        z_res[s,] <- mat_z[genes, ]
    } else {
        warning(sprintf("No genes found for signature %s", s))
    }
}

if(SCORING_METHOD == "both") {
    cor_methods <- cor(as.vector(gsva_res), as.vector(z_res), use="complete.obs")
    
    add_stat_log("Scoring Method Comparison", sprintf(
        "GSVA vs Z-Score correlation: r = %.3f\n  Interpretation: %s\n  Decision: Using Z-Score (more robust for small gene sets and small N)",
        cor_methods,
        if(cor_methods > 0.8) "High agreement" else if(cor_methods > 0.6) "Moderate agreement" else "Low agreement - interpret with caution"
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
        theme_publication(base_size=12)
    
    ggsave(paste0(opt$out, "_Method_Agreement_mqc.png"), p_agreement, width=8, height=6)
    
    final_scores <- z_res
    scoring_label <- "Z-Score (Robust)"
} else if(SCORING_METHOD == "zscore") {
    final_scores <- z_res
    scoring_label <- "Z-Score"
    add_stat_log("Scoring Method", "Using Z-Score method (user-specified)")
} else {
    final_scores <- gsva_res
    scoring_label <- "GSVA"
    add_stat_log("Scoring Method", "Using GSVA method (user-specified)")
}

cat(sprintf("→ Using %s for downstream analysis\n", scoring_label))

# ==============================================================================
# 4. PLASTICITY (SHANNON ENTROPY WITH DUAL TESTING - BOTH TESTS BY DEFAULT)
# ==============================================================================
cat("LOG [4/11]: Plasticity Analysis (running BOTH ANOVA and Kruskal-Wallis)...\n")

calc_entropy <- function(s) {
    s <- pmin(pmax(s, -10), 10)
    p <- exp(s)/sum(exp(s))
    -sum(p*log(p + 1e-10), na.rm=TRUE)
}

meta$Plasticity <- apply(final_scores, 2, calc_entropy)

n_total <- nrow(meta)

# Normality test
shapiro_test <- shapiro.test(meta$Plasticity)
shapiro_p <- shapiro_test$p.value
shapiro_w <- shapiro_test$statistic

# ALWAYS run BOTH tests for comparison
cat("  → Running ANOVA (parametric)...\n")
plast_aov <- aov(Plasticity ~ Classification, data=meta)
plast_summary <- summary(plast_aov)[[1]]
anova_p <- plast_summary[["Pr(>F)"]][1]
anova_f <- plast_summary[["F value"]][1]
anova_df1 <- plast_summary[["Df"]][1]
anova_df2 <- plast_summary[["Df"]][2]

cat("  → Running Kruskal-Wallis (non-parametric)...\n")
kw_test <- kruskal.test(Plasticity ~ Classification, data=meta)
kw_p <- kw_test$p.value
kw_h <- kw_test$statistic
kw_df <- kw_test$parameter

# Create comparison table
plasticity_tests <- data.frame(
    Test = c("Shapiro-Wilk (Normality)", "ANOVA (Parametric)", "Kruskal-Wallis (Non-parametric)"),
    Statistic = c(
        sprintf("W = %.3f", shapiro_w),
        sprintf("F(%d,%d) = %.3f", anova_df1, anova_df2, anova_f),
        sprintf("H(%d) = %.3f", kw_df, kw_h)
    ),
    P_value = c(safe_format(shapiro_p), safe_format(anova_p), safe_format(kw_p)),
    Significance = c(interpret_p(shapiro_p), interpret_p(anova_p), interpret_p(kw_p)),
    Interpretation = c(
        if(shapiro_p >= 0.05) "Data is normally distributed" else "Data is non-normal",
        if(anova_p < 0.05) "Groups differ significantly" else "No significant differences",
        if(kw_p < 0.05) "Groups differ significantly" else "No significant differences"
    ),
    stringsAsFactors = FALSE
)

# Log both results
add_stat_log("Plasticity: Normality Test", sprintf(
    "Shapiro-Wilk W = %.3f, P = %s %s\n  Interpretation: Data is %s",
    shapiro_w, safe_format(shapiro_p), interpret_p(shapiro_p),
    if(shapiro_p >= 0.05) "normally distributed" else "non-normal"
))

add_stat_log("Plasticity: ANOVA (Parametric)", sprintf(
    "One-way ANOVA F(%d,%d) = %.3f\n  P = %s %s\n  Interpretation: %s",
    anova_df1, anova_df2, anova_f,
    safe_format(anova_p), interpret_p(anova_p),
    if(anova_p < 0.05) "Significant differences in plasticity" else "No significant differences"
))

add_stat_log("Plasticity: Kruskal-Wallis (Non-parametric)", sprintf(
    "Kruskal-Wallis H(%d) = %.3f\n  P = %s %s\n  Interpretation: %s",
    kw_df, kw_h,
    safe_format(kw_p), interpret_p(kw_p),
    if(kw_p < 0.05) "Significant differences in plasticity" else "No significant differences"
))

# Determine which test to use for plot subtitle
if(n_total < 30 || shapiro_p < 0.05) {
    primary_test <- "Kruskal-Wallis"
    primary_p <- kw_p
    cat(sprintf("  → Recommendation: Use Kruskal-Wallis (N=%d, normality P=%s)\n", n_total, safe_format(shapiro_p)))
} else {
    primary_test <- "ANOVA"
    primary_p <- anova_p
    cat(sprintf("  → Recommendation: Use ANOVA (N=%d, normality P=%s)\n", n_total, safe_format(shapiro_p)))
}

plast_stats <- meta %>%
    group_by(Classification) %>%
    summarise(
        N = n(),
        Mean = mean(Plasticity),
        SD = sd(Plasticity),
        Median = median(Plasticity),
        IQR = IQR(Plasticity),
        .groups="drop"
    )

add_stat_log("Plasticity: Descriptive Statistics", sprintf(
    "Group statistics:\n  %s",
    paste(capture.output(print(plast_stats, row.names=FALSE)), collapse="\n  ")
))

# Export plasticity test comparison
write.csv(plasticity_tests, paste0(opt$out, "_Plasticity_Tests_Comparison.csv"), row.names=FALSE)

p_plast <- ggplot(meta, aes(x=Classification, y=Plasticity, fill=Classification)) +
    geom_violin(alpha=0.3, trim=FALSE) +
    geom_boxplot(width=0.2, fill="white", outlier.shape=NA) +
    geom_jitter(width=0.1, size=3, alpha=0.7) +
    stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="white") +
    scale_fill_manual(values=GROUP_COLORS) +
    labs(title="Cellular Plasticity (Shannon Entropy)",
         subtitle=sprintf("Recommended: %s P=%s %s | Both tests shown in report", 
                         primary_test, safe_format(primary_p), interpret_p(primary_p)),
         caption="Diamond = mean, box = median ± IQR, dots = individual samples") +
    theme_publication() +
    theme(legend.position="none")

ggsave(paste0(opt$out, "_Plasticity_mqc.png"), p_plast, width=7, height=6)

# ==============================================================================
# 5. DUAL-WEIGHTING DIFFERENTIAL ANALYSIS (AUTOMATIC COMPARISON!)
# ==============================================================================
cat("LOG [5/11]: Dual-Weighting Differential Analysis (weighted vs unweighted)...\n")
cat("  → Running BOTH analyses automatically for comparison...\n\n")

# Function to run limma analysis with optional weights
run_limma_analysis <- function(scores, metadata, use_weights=TRUE, label="") {
    cat(sprintf("  [%s] Running analysis...\n", label))
    
    design <- model.matrix(~0 + metadata$Classification)
    colnames(design) <- levels(metadata$Classification)
    
    if(use_weights) {
        aw <- arrayWeights(scores, design)
        weight_df <- data.frame(
            Sample = names(aw),
            Weight = round(aw, 3),
            Group = as.character(metadata[names(aw), "Classification"])
        )
        outlier_threshold <- mean(aw) - 2*sd(aw)
        outliers <- weight_df$Sample[weight_df$Weight < outlier_threshold]
        
        cat(sprintf("    arrayWeights: mean=%.3f, range=%.3f-%.3f\n", 
                   mean(aw), min(aw), max(aw)))
        if(length(outliers) > 0) {
            cat(sprintf("    ⚠ Potential outliers: %s\n", paste(outliers, collapse=", ")))
        }
    } else {
        aw <- rep(1, ncol(scores))
        names(aw) <- colnames(scores)
        weight_df <- data.frame(
            Sample = names(aw),
            Weight = round(aw, 3),
            Group = as.character(metadata[names(aw), "Classification"])
        )
        cat("    Equal weights: all samples = 1.000\n")
    }
    
    # Standard pairwise contrasts
    levs <- levels(metadata$Classification)
    contrast_formulas <- c()
    contrast_names <- c()
    
    if(length(levs) >= 2) {
        contrast_formulas <- c(contrast_formulas, sprintf("%s - %s", levs[2], levs[1]))
        contrast_names <- c(contrast_names, sprintf("%s_vs_%s", 
                                                    gsub("[^A-Za-z0-9]", "", levs[2]), 
                                                    gsub("[^A-Za-z0-9]", "", levs[1])))
    }
    if(length(levs) >= 3) {
        contrast_formulas <- c(contrast_formulas,
                              sprintf("%s - %s", levs[3], levs[1]),
                              sprintf("%s - %s", levs[3], levs[2]))
        contrast_names <- c(contrast_names,
                           sprintf("%s_vs_%s", gsub("[^A-Za-z0-9]", "", levs[3]), gsub("[^A-Za-z0-9]", "", levs[1])),
                           sprintf("%s_vs_%s", gsub("[^A-Za-z0-9]", "", levs[3]), gsub("[^A-Za-z0-9]", "", levs[2])))
    }
    
    cont.matrix <- makeContrasts(contrasts=contrast_formulas, levels=design)
    colnames(cont.matrix) <- contrast_names
    
    # NEW: Polynomial contrasts for trajectory testing
    # Test for LINEAR and QUADRATIC trends across ordered groups
    design_poly <- model.matrix(~ poly(as.numeric(metadata$Classification), 2))
    colnames(design_poly) <- c("Intercept", "Linear", "Quadratic")
    
    # Fit models
    fit <- lmFit(scores, design, weights=aw)
    fit <- contrasts.fit(fit, cont.matrix)
    fit <- eBayes(fit)
    
    fit_poly <- lmFit(scores, design_poly, weights=aw)
    fit_poly <- eBayes(fit_poly)
    
    return(list(
        fit = fit,
        fit_poly = fit_poly,
        weights = weight_df,
        contrasts = contrast_names,
        levs = levs
    ))
}

# Run BOTH analyses
cat("\n╔═══════════════════════════════════════════════╗\n")
cat("║   ANALYSIS 1: WITH arrayWeights (Standard)   ║\n")
cat("╚═══════════════════════════════════════════════╝\n")
res_weighted <- run_limma_analysis(final_scores, meta, use_weights=TRUE, label="WEIGHTED")

cat("\n╔═══════════════════════════════════════════════╗\n")
cat("║  ANALYSIS 2: WITHOUT arrayWeights (Unbiased)  ║\n")
cat("╚═══════════════════════════════════════════════╝\n")
res_unweighted <- run_limma_analysis(final_scores, meta, use_weights=FALSE, label="UNWEIGHTED")

# Store for later comparison
all_coefs <- res_weighted$contrasts
levs <- res_weighted$levs

# ==============================================================================
# 6. TRAJECTORY TREND TESTING (NEW: Jonckheere-Terpstra!)
# ==============================================================================
cat("\nLOG [6/11]: Trajectory Trend Testing (Jonckheere-Terpstra)...\n")
cat("  → Testing for monotonic trends across Culture → Primary → Recurrent\n\n")

stage_numeric <- as.numeric(meta$Classification)

# Initialize results dataframe
trajectory_results <- data.frame(
    Signature = rownames(final_scores),
    stringsAsFactors = FALSE
)

# Run Jonckheere-Terpstra test for each signature - FOR BOTH Z-SCORE AND GSVA
for(method_name in c("ZScore", "GSVA")) {
    cat(sprintf("  → Running JT tests for %s...\n", method_name))
    curr_scores <- if(method_name == "ZScore") z_res else gsva_res
    
    for(i in 1:nrow(curr_scores)) {
        sig_name <- rownames(curr_scores)[i]
        scores <- curr_scores[i, ]
        
        # JT test with permutation-based p-value
        jt <- tryCatch({
            jonckheere.test(scores, stage_numeric, nperm=1000)
        }, error = function(e) {
            list(statistic = NA, p.value = NA)
        })
        
        trajectory_results[[paste0(method_name, "_JT_Statistic")]][i] <- jt$statistic
        trajectory_results[[paste0(method_name, "_JT_P")]][i] <- jt$p.value
        
        # Determine direction
        n <- length(scores)
        g_sizes <- table(stage_numeric)
        expected_mean <- (n^2 - sum(g_sizes^2)) / 4
        
        trajectory_results[[paste0(method_name, "_JT_Direction")]][i] <- ifelse(
            is.na(jt$statistic), "Unknown",
            ifelse(jt$statistic > expected_mean, "Increasing", "Decreasing")
        )
        
        trajectory_results[[paste0(method_name, "_JT_Sig")]][i] <- interpret_p(jt$p.value)
    }
}

# Add polynomial contrast results from BOTH weighting methods
for(method_name in c("Weighted", "Unweighted")) {
    fit_poly <- if(method_name == "Weighted") res_weighted$fit_poly else res_unweighted$fit_poly
    
    tt_linear <- topTable(fit_poly, coef="Linear", number=Inf, sort.by="none")
    tt_quad <- topTable(fit_poly, coef="Quadratic", number=Inf, sort.by="none")
    
    trajectory_results[[paste0(method_name, "_Linear_P")]] <- tt_linear[trajectory_results$Signature, "adj.P.Val"]
    trajectory_results[[paste0(method_name, "_Linear_Coef")]] <- tt_linear[trajectory_results$Signature, "logFC"]
    trajectory_results[[paste0(method_name, "_Quad_P")]] <- tt_quad[trajectory_results$Signature, "adj.P.Val"]
    trajectory_results[[paste0(method_name, "_Quad_Coef")]] <- tt_quad[trajectory_results$Signature, "logFC"]
}

# Classify trajectory patterns (using Z-Score as primary)
trajectory_results$Pattern <- sapply(1:nrow(trajectory_results), function(i) {
    lin_p_w <- trajectory_results$Weighted_Linear_P[i]
    quad_p_w <- trajectory_results$Weighted_Quad_P[i]
    jt_p <- trajectory_results$ZScore_JT_P[i]
    
    if(is.na(jt_p)) return("No trend")
    
    if(jt_p < 0.05 && lin_p_w < 0.05 && quad_p_w >= 0.05) {
        return("Linear (monotonic)")
    } else if(quad_p_w < 0.05) {
        return("Quadratic (spike/dip)")
    } else if(jt_p < 0.1) {
        return("Weak trend")
    } else {
        return("No trend")
    }
})

# Determine consensus
trajectory_results$Consensus_Trend <- sapply(1:nrow(trajectory_results), function(i) {
    z_sig <- !is.na(trajectory_results$ZScore_JT_P[i]) && trajectory_results$ZScore_JT_P[i] < 0.05
    g_sig <- !is.na(trajectory_results$GSVA_JT_P[i]) && trajectory_results$GSVA_JT_P[i] < 0.05
    
    if(z_sig && g_sig) {
        return("ROBUST (Both)")
    } else if (z_sig) {
        return("Z-Score Only")
    } else if (g_sig) {
        return("GSVA Only")
    } else {
        return("None")
    }
})

cat("\n📊 Trajectory Pattern Summary:\n")
print(table(trajectory_results$Pattern))
cat("\n📊 Consensus Summary:\n")
print(table(trajectory_results$Consensus_Trend))

add_stat_log("Trajectory Testing", sprintf(
    "Jonckheere-Terpstra trend analysis:\n  %s",
    paste(capture.output(print(table(trajectory_results$Pattern))), collapse="\n  ")
))

# Export trajectory results
write.csv(trajectory_results, paste0(opt$out, "_Trajectory_Tests.csv"), row.names=FALSE)

# ==============================================================================
# 7. MULTI-FDR ANALYSIS + WEIGHTING COMPARISON
# ==============================================================================
cat("\nLOG [7/11]: Multi-FDR Analysis with Weighting Comparison...\n")

# Compute signature correlations
sig_cor <- cor(t(final_scores), method="pearson")

cor_pairs <- which(abs(sig_cor) > 0.6 & upper.tri(sig_cor), arr.ind=TRUE)
if(nrow(cor_pairs) > 0) {
    cor_summary <- data.frame(
        Sig1 = rownames(sig_cor)[cor_pairs[,1]],
        Sig2 = rownames(sig_cor)[cor_pairs[,2]],
        Correlation = sig_cor[cor_pairs]
    )
    cor_summary <- cor_summary[order(abs(cor_summary$Correlation), decreasing=TRUE), ]
} else {
    cor_summary <- data.frame()
}

# Initialize ENHANCED LLM summary
llm_summary <- paste0(
    "═══════════════════════════════════════════════════════════════════\n",
    "GBM LITT THERAPY SUBTYPE EVOLUTION ANALYSIS (OMNIBUS EDITION)\n",
    "═══════════════════════════════════════════════════════════════════\n\n",
    "EXPERIMENTAL DESIGN:\n",
    sprintf("  Trajectory: %s\n", paste(levs, collapse=" → ")),
    sprintf("  Samples: %s\n", paste(names(group_sizes), "=", group_sizes, collapse=", ")),
    sprintf("  Total signatures tested: %d\n", length(sigs)),
    sprintf("  Scoring methods: GSVA + Z-Score (dual approach)\n"),
    sprintf("  Research Question: Does Litt therapy induce subtype evolution?\n\n"),
    
    "═══════════════════════════════════════════════════════════════════\n",
    "1. GLOBAL TRANSCRIPTOMIC STRUCTURE (PCA)\n",
    "═══════════════════════════════════════════════════════════════════\n",
    sprintf("  PC1 variance explained: %.1f%%\n", var_pc[1]),
    sprintf("  PC2 variance explained: %.1f%%\n", var_pc[2]),
    sprintf("  Cumulative (PC1+PC2): %.1f%%\n\n", cum_var[2]),
    
    "  PERMANOVA (Group Separation Test):\n",
    sprintf("    F-statistic: %.3f\n", perm_f),
    sprintf("    R² (variance explained by groups): %.3f (%.1f%%)\n", perm_r2, perm_r2*100),
    sprintf("    P-value: %s %s\n", safe_format(perm_p), interpret_p(perm_p)),
    sprintf("    Interpretation: Stages %s transcriptionally distinct\n\n",
            if(!is.na(perm_p) && perm_p < 0.05) "ARE" else "are NOT"),
    
    "  Top Gene Drivers by Principal Component:\n"
)

for(pc in names(pc_drivers_list)) {
    llm_summary <- paste0(llm_summary, sprintf("    %s: %s\n", pc, paste(pc_drivers_list[[pc]], collapse=", ")))
}

llm_summary <- paste0(llm_summary, "\n",
    "═══════════════════════════════════════════════════════════════════\n",
    "2. CELLULAR PLASTICITY (Shannon Entropy)\n",
    "═══════════════════════════════════════════════════════════════════\n",
    sprintf("  Normality test (Shapiro-Wilk): W=%.3f, P=%s %s\n",
            shapiro_w, safe_format(shapiro_p), interpret_p(shapiro_p)),
    sprintf("  ANOVA: F(%.0f,%.0f)=%.3f, P=%s %s\n",
            anova_df1, anova_df2, anova_f, safe_format(anova_p), interpret_p(anova_p)),
    sprintf("  Kruskal-Wallis: H(%.0f)=%.3f, P=%s %s\n\n",
            kw_df, kw_h, safe_format(kw_p), interpret_p(kw_p)),
    
    "  Group Statistics:\n"
)

for(i in 1:nrow(plast_stats)) {
    llm_summary <- paste0(llm_summary,
        sprintf("    %s: Mean=%.3f±%.3f, Median=%.3f, IQR=%.3f (n=%d)\n",
                plast_stats$Classification[i], plast_stats$Mean[i], plast_stats$SD[i],
                plast_stats$Median[i], plast_stats$IQR[i], plast_stats$N[i]))
}

llm_summary <- paste0(llm_summary, "\n",
    "═══════════════════════════════════════════════════════════════════\n",
    "3. SIGNATURE CO-EVOLUTION (Correlation Analysis)\n",
    "═══════════════════════════════════════════════════════════════════\n"
)

if(nrow(cor_summary) > 0) {
    llm_summary <- paste0(llm_summary, sprintf("  High Correlations (|r| > 0.6): %d pairs detected\n\n", nrow(cor_summary)))
    for(i in 1:min(10, nrow(cor_summary))) {
        llm_summary <- paste0(llm_summary,
            sprintf("    %s ↔ %s: r=%.3f\n",
                    cor_summary$Sig1[i], cor_summary$Sig2[i], cor_summary$Correlation[i]))
    }
    if(nrow(cor_summary) > 10) {
        llm_summary <- paste0(llm_summary, sprintf("    ... and %d more (see correlation matrix)\n", nrow(cor_summary) - 10))
    }
} else {
    llm_summary <- paste0(llm_summary, "  No strong correlations (|r| > 0.6) detected between signatures\n")
}

llm_summary <- paste0(llm_summary, "\n",
    "═══════════════════════════════════════════════════════════════════\n",
    "4. TRAJECTORY TRENDS (Dual Method Consensus)\n",
    "═══════════════════════════════════════════════════════════════════\n",
    "Testing for monotonic trends using Jonckheere-Terpstra test\n",
    "Methods: Z-Score + GSVA for consensus confirmation\n\n"
)

sig_trends <- trajectory_results[trajectory_results$Consensus_Trend != "None", ]
if(nrow(sig_trends) > 0) {
    sig_trends <- sig_trends[order(match(sig_trends$Consensus_Trend, c("ROBUST (Both)", "Z-Score Only", "GSVA Only"))), ]
    for(i in 1:nrow(sig_trends)) {
        sig <- rownames(sig_trends)[i]
        llm_summary <- paste0(llm_summary, sprintf(
            "  %s [%s]:\n    Z-Score: P=%s (Stat=%.1f, Direction=%s)\n    GSVA: P=%s (Stat=%.1f, Direction=%s)\n    Pattern: %s\n\n",
            sig, sig_trends$Consensus_Trend[i],
            safe_format(sig_trends$ZScore_JT_P[i]), sig_trends$ZScore_JT_Statistic[i], sig_trends$ZScore_JT_Direction[i],
            safe_format(sig_trends$GSVA_JT_P[i]), sig_trends$GSVA_JT_Statistic[i], sig_trends$GSVA_JT_Direction[i],
            sig_trends$Pattern[i]
        ))
    }
} else {
    llm_summary <- paste0(llm_summary, "  No significant monotonic trends detected in either method.\n")
}

llm_summary <- paste0(llm_summary, "\n═══════════════════════════════════════════════════════════════════\n",
    "DIFFERENTIAL ANALYSIS: WEIGHTED vs UNWEIGHTED COMPARISON\n",
    "═══════════════════════════════════════════════════════════════════\n",
    "(Weighted = arrayWeights ON; Unweighted = all samples equal)\n\n"
)

# Run multi-FDR for BOTH weighting methods
for(method_name in c("Weighted", "Unweighted")) {
    cat(sprintf("\n  → Processing %s analysis...\n", method_name))
    
    fit_current <- if(method_name == "Weighted") res_weighted$fit else res_unweighted$fit
    
    llm_summary <- paste0(llm_summary, 
                         sprintf("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"),
                         sprintf("%s ANALYSIS\n", toupper(method_name)),
                         sprintf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"))
    
    for(thresh in FDR_THRESHOLDS) {
        llm_summary <- paste0(llm_summary, sprintf("\nFDR < %.3f:\n", thresh))
        
        for(coef in all_coefs) {
            tt <- topTable(fit_current, coef=coef, number=Inf, adjust.method="BH")
            sig_hits <- rownames(tt)[tt$adj.P.Val < thresh]
            
            llm_summary <- paste0(llm_summary, sprintf("\n  %s:\n", coef))
            if(length(sig_hits) > 0) {
                for(sig in sig_hits) {
                    direction <- if(tt[sig, "logFC"] > 0) "UP" else "DOWN"
                    llm_summary <- paste0(llm_summary,
                        sprintf("    • %s: %s (logFC=%.3f, FDR=%s)\n",
                                sig, direction, tt[sig, "logFC"], safe_format(tt[sig, "adj.P.Val"])))
                }
            } else {
                llm_summary <- paste0(llm_summary, "    • No significant changes\n")
            }
        }
    }
}

llm_summary <- paste0(llm_summary, "\n\n",
    "═══════════════════════════════════════════════════════════════════\n",
    "INTERPRETATION GUIDE FOR LLM:\n",
    "═══════════════════════════════════════════════════════════════════\n",
    "Consensus Strength:\n",
    "  - ROBUST (Both): Significant in both Z-Score AND GSVA → High confidence\n",
    "  - Method-specific: Significant in only one method → Moderate confidence\n",
    "  - Neither: Not significant → No evidence of trend\n\n",
    
    "Compare WEIGHTED vs UNWEIGHTED results:\n",
    "  1. If a signature is significant ONLY in unweighted → arrayWeights masked it\n",
    "  2. If significant in BOTH → robust finding (high confidence)\n",
    "  3. Linear patterns = steady evolution (e.g., EMT progression)\n",
    "  4. Quadratic patterns = transient state (e.g., therapy shock response)\n\n",
    
    "Focus on signatures with:\n",
    "  - Jonckheere-Terpstra P < 0.05 (true evolutionary trend)\n",
    "  - Consensus across BOTH scoring methods (robust)\n",
    "  - Appearing in BOTH weighting methods (robust)\n",
    "  - Clinical relevance to GBM recurrence mechanisms\n\n",
    
    "Significance Codes:\n",
    "  *** P < 0.001 (highly significant)\n",
    "  **  P < 0.01  (very significant)\n",
    "  *   P < 0.05  (significant)\n",
    "  .   P < 0.10  (trend)\n",
    "  ns  P ≥ 0.10  (not significant)\n",
    "═══════════════════════════════════════════════════════════════════\n"
)

llm_file <- paste0(opt$out, "_llm_summary.txt")
writeLines(llm_summary, llm_file)
cat(sprintf("\n✓ LLM interpretation prompt: %s\n", basename(llm_file)))

# ==============================================================================
# 8. UNIFIED TRAJECTORY PLOT (FIXED - ALL SIGNATURES SHOW LINES!)
# ==============================================================================
cat("\nLOG [8/11]: Trajectory Analysis (with fixed line rendering)...\n")

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

# *** THE FIX: na.rm=TRUE ***
traj_summary <- traj_data %>%
    group_by(Signature, Class, Stage) %>%
    summarise(Mean = mean(Score, na.rm=TRUE), SE = sd(Score, na.rm=TRUE)/sqrt(n()), .groups="drop")

trend_order <- trajectory_results %>%
    arrange(desc(ZScore_JT_Statistic)) %>%
    pull(Signature)

traj_data$Signature <- factor(traj_data$Signature, levels=trend_order)
traj_summary$Signature <- factor(traj_summary$Signature, levels=trend_order)

p_traj <- ggplot(traj_data, aes(x=Stage, y=Score)) +
    geom_ribbon(data=traj_summary, 
                aes(x=Stage, ymin=Mean-SE, ymax=Mean+SE, fill=Signature, group=Signature),
                inherit.aes=FALSE, alpha=0.2) +
    geom_line(data=traj_summary, 
              aes(x=Stage, y=Mean, color=Signature, group=Signature),
              inherit.aes=FALSE, linewidth=1.2, alpha=0.9) +
    geom_point(aes(fill=Class, shape=Class), size=3, alpha=0.6) +
    facet_wrap(~Signature, scales="free_y", ncol=3) +
    scale_x_continuous(breaks=1:length(levs), labels=levs) +
    scale_fill_manual(values=GROUP_COLORS, name="Stage") +
    scale_color_brewer(palette="Set1", guide="none") +
    scale_shape_manual(values=GROUP_SHAPES, name="Stage") +
    labs(title="Signature Trajectories Across Litt Therapy Evolution",
         subtitle="Ordered by JT trend statistic (top = strongest increasing) | Lines = group means ± SE | FIXED: All lines now render",
         x="Stage", y="Z-Score Expression") +
    theme_publication(base_size=11) +
    theme(legend.position = "bottom")

ggsave(paste0(opt$out, "_Unified_Trajectories_mqc.png"), p_traj, width=14, height=12)

# ==============================================================================
# 9. HEATMAPS & CORRELATION
# ==============================================================================
cat("LOG [9/11]: Heatmaps & Co-evolution Analysis...\n")

ha <- HeatmapAnnotation(
    Stage = meta$Classification,
    Plasticity = meta$Plasticity,
    col = list(
        Stage = GROUP_COLORS,
        Plasticity = colorRamp2(c(min(meta$Plasticity), max(meta$Plasticity)), c("white", "purple"))
    )
)

png(paste0(opt$out, "_Heatmap_mqc.png"), width=12, height=8, units="in", res=300)
Heatmap(final_scores, name="Z-Score", top_annotation=ha,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        cluster_columns=FALSE, column_order=order(meta$Classification),
        row_names_gp = gpar(fontsize=9))
dev.off()

png(paste0(opt$out, "_Sig_Correlation_mqc.png"), width=8, height=7, units="in", res=300)
Heatmap(sig_cor, name="Pearson\nCorr",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        column_title = "Signature Co-evolution",
        row_names_gp = gpar(fontsize=10),
        column_names_gp = gpar(fontsize=10))
dev.off()

# ==============================================================================
# 10. EXPORT DATA FILES
# ==============================================================================
cat("LOG [10/11]: Exporting Data Files...\n")

write.csv(t(final_scores), paste0(opt$out, "_Scores.csv"))
write.csv(t(z_res), paste0(opt$out, "_ZScores.csv"))
write.csv(t(gsva_res), paste0(opt$out, "_GSVA_Scores.csv"))
write.csv(meta, paste0(opt$out, "_Metadata.csv"))
write.csv(res_weighted$weights, paste0(opt$out, "_Weights_Weighted.csv"), row.names=FALSE)
write.csv(res_unweighted$weights, paste0(opt$out, "_Weights_Unweighted.csv"), row.names=FALSE)
write.csv(plast_stats, paste0(opt$out, "_Plasticity_Statistics.csv"), row.names=FALSE)
write.csv(trajectory_results, paste0(opt$out, "_Trajectory_Tests.csv"), row.names=FALSE)
if(nrow(cor_summary) > 0) {
    write.csv(cor_summary, paste0(opt$out, "_Signature_Correlations.csv"), row.names=FALSE)
}

# Export MultiQC-compatible summary
summary_file <- paste0(dirname(opt$out), "/trajectory_summary_mqc.tsv")
cat("# id: 'trajectory_tests'\n# section_name: 'Trajectory Trend Tests'\n# plot_type: 'table'\n", 
    file=summary_file)
write.table(trajectory_results, file=summary_file, sep="\t", quote=FALSE, row.names=FALSE, append=TRUE)

# ==============================================================================
# 11. HTML REPORT (WITH EMBEDDED LLM SUMMARY)
# ==============================================================================
cat("LOG [11/11]: Generating HTML Report...\n")

summary_html <- paste0(dirname(opt$out), "/analysis_summary_mqc.html")
sink(summary_html)

cat("
<div style='background-color:#f8f9fa; padding:20px; border:2px solid #667eea; border-radius:8px;'>
    <h2 style='color:#667eea;'>🧬 Litt Therapy Subtype Evolution Analysis (v16.2 Omnibus)</h2>
    
    <div style='background:#fff3cd; border-left:4px solid #ffc107; padding:10px; margin:15px 0;'>
    <strong>🆕 v16.2 Features:</strong><br>
    ✓ Automatic dual-weighting comparison (arrayWeights ON vs OFF)<br>
    ✓ Jonckheere-Terpstra trajectory testing (monotonic trends)<br>
    ✓ Polynomial contrasts (linear/quadratic patterns)<br>
    ✓ Both ANOVA and Kruskal-Wallis plasticity tests<br>
    ✓ Dual scoring (GSVA + Z-Score) with consensus reporting<br>
    ✓ Enhanced LLM summary with PCA, plasticity, and correlations<br>
    ✓ Fixed trajectory plot rendering (all lines now show)
    </div>
    
    <h3>Dataset Overview</h3>
    <pre>", paste(names(group_sizes), ":", group_sizes, collapse="\n"), "</pre>
    
    <h3>Global Structure</h3>
    <ul>
        <li><strong>PERMANOVA:</strong> R²=", sprintf("%.3f", perm_r2), " (", sprintf("%.1f", perm_r2*100), 
            "% variance explained), P=", safe_format(perm_p), " ", interpret_p(perm_p), "</li>
        <li><strong>PCA:</strong> PC1=", var_pc[1], "%, PC2=", var_pc[2], "%, Total=", cum_var[2], "%</li>
    </ul>
    
    <h3>Scoring Method Agreement</h3>
    <ul>
        <li><strong>GSVA vs Z-Score correlation:</strong> r=", sprintf("%.3f", cor_methods), 
            " (", if(cor_methods > 0.8) "High agreement" else if(cor_methods > 0.6) "Moderate agreement" else "Low agreement", ")</li>
    </ul>
    
    <h3>Plasticity Analysis (Dual Testing)</h3>
    <p style='color:#666; font-size:13px;'>
    <strong>Why run both tests?</strong> ANOVA is more powerful when data is normal, but Kruskal-Wallis is more robust to outliers and non-normality. 
    For small sample sizes (N<30), Kruskal-Wallis is often preferred. When results agree, confidence is higher.
    </p>
    <table style='border-collapse:collapse; width:100%; margin:15px 0; border:1px solid #ddd;'>
        <tr style='background:#667eea; color:white;'>
            <th style='padding:8px; border:1px solid #555;'>Test</th>
            <th style='padding:8px; border:1px solid #555;'>Statistic</th>
            <th style='padding:8px; border:1px solid #555;'>P-value</th>
            <th style='padding:8px; border:1px solid #555;'>Sig</th>
            <th style='padding:8px; border:1px solid #555;'>Interpretation</th>
        </tr>
")

for(i in 1:nrow(plasticity_tests)) {
    row_color <- if(i == 1) "#f0f0f0" else "white"
    cat(sprintf("        <tr style='background:%s;'>
            <td style='padding:8px; border:1px solid #ddd;'><strong>%s</strong></td>
            <td style='padding:8px; border:1px solid #ddd;'>%s</td>
            <td style='padding:8px; border:1px solid #ddd;'>%s</td>
            <td style='padding:8px; border:1px solid #ddd;'>%s</td>
            <td style='padding:8px; border:1px solid #ddd;'>%s</td></tr>\n",
            row_color,
            plasticity_tests$Test[i],
            plasticity_tests$Statistic[i],
            plasticity_tests$P_value[i],
            plasticity_tests$Significance[i],
            plasticity_tests$Interpretation[i]))
}

cat("    </table>
    <p style='font-size:12px; color:#666; font-style:italic;'>
    <strong>Recommendation:</strong> ", 
    if(n_total < 30 || shapiro_p < 0.05) {
        sprintf("Use <strong>Kruskal-Wallis</strong> (N=%d, data is %s)", 
                n_total, if(shapiro_p < 0.05) "non-normal" else "small sample")
    } else {
        sprintf("Use <strong>ANOVA</strong> (N=%d, data is normally distributed)", n_total)
    },
    ". When both tests agree (both significant or both non-significant), confidence is highest.
    </p>
    
    <h3>Plasticity Descriptive Statistics by Group</h3>
    <table style='border-collapse: collapse; width:100%; margin:15px 0;'>
        <tr style='background:#667eea; color:white;'>
            <th style='padding:8px; border:1px solid #555;'>Group</th>
            <th style='padding:8px; border:1px solid #555;'>N</th>
            <th style='padding:8px; border:1px solid #555;'>Mean</th>
            <th style='padding:8px; border:1px solid #555;'>SD</th>
            <th style='padding:8px; border:1px solid #555;'>Median</th>
            <th style='padding:8px; border:1px solid #555;'>IQR</th>
        </tr>
")

for(i in 1:nrow(plast_stats)) {
    cat(sprintf("        <tr><td style='padding:8px; border:1px solid #ccc;'>%s</td>
            <td style='padding:8px; border:1px solid #ccc;'>%d</td>
            <td style='padding:8px; border:1px solid #ccc;'>%.3f</td>
            <td style='padding:8px; border:1px solid #ccc;'>%.3f</td>
            <td style='padding:8px; border:1px solid #ccc;'>%.3f</td>
            <td style='padding:8px; border:1px solid #ccc;'>%.3f</td></tr>\n",
            plast_stats$Classification[i], plast_stats$N[i], plast_stats$Mean[i],
            plast_stats$SD[i], plast_stats$Median[i], plast_stats$IQR[i]))
}

cat("    </table>
    
    <h3>Trajectory Consensus Summary</h3>
    <table style='border-collapse:collapse; width:100%;'>
        <tr style='background:#667eea; color:white;'>
            <th style='padding:8px;'>Consensus</th>
            <th style='padding:8px;'>Count</th>
        </tr>
")

consensus_counts <- table(trajectory_results$Consensus_Trend)
for(cons in names(consensus_counts)) {
    cat(sprintf("<tr><td style='padding:8px;'>%s</td><td style='padding:8px;'>%d</td></tr>\n",
               cons, consensus_counts[cons]))
}

cat("    </table>
    
    <h3>Top Trajectory Signatures (Consensus Trends)</h3>
")

if(nrow(sig_trends) > 0) {
    cat("<table style='border-collapse:collapse; width:100%;'>
        <tr style='background:#667eea; color:white;'>
            <th style='padding:8px;'>Signature</th>
            <th style='padding:8px;'>Consensus</th>
            <th style='padding:8px;'>Z-Score P</th>
            <th style='padding:8px;'>GSVA P</th>
            <th style='padding:8px;'>Pattern</th>
        </tr>")
    
    for(i in 1:nrow(sig_trends)) {
        cat(sprintf("<tr><td style='padding:8px;'><strong>%s</strong></td>
                    <td style='padding:8px;'>%s</td>
                    <td style='padding:8px;'>%s</td>
                    <td style='padding:8px;'>%s</td>
                    <td style='padding:8px;'>%s</td></tr>\n",
                   rownames(sig_trends)[i], sig_trends$Consensus_Trend[i],
                   safe_format(sig_trends$ZScore_JT_P[i]),
                   safe_format(sig_trends$GSVA_JT_P[i]),
                   sig_trends$Pattern[i]))
    }
    cat("</table>")
} else {
    cat("<p>No significant monotonic trends detected.</p>")
}

cat("
</div>

<div style='background-color:#e7f3ff; padding:20px; border:2px solid #667eea; border-radius:8px; margin-top:20px;'>
    <h2 style='color:#667eea; margin-top:0;'>📄 Enhanced LLM Interpretation Summary</h2>
    <p style='color:#666; font-size:13px;'>
    <strong>Note:</strong> This summary is also available as a standalone file: <code>", basename(llm_file), "</code><br>
    Copy the text below into ChatGPT/Claude/Gemini for biological interpretation.<br>
    <strong>NEW:</strong> Now includes PCA components, plasticity statistics, correlation analysis, and dual-method consensus.
    </p>
    
    <div style='background:#fff; padding:15px; border:1px solid #ddd; border-radius:5px; margin-top:15px;'>
    <pre style='font-family:monospace; font-size:11px; white-space:pre-wrap; margin:0;'>",
    llm_summary,
    "</pre>
    </div>
</div>

<div style='background:#fff; padding:15px; border:1px solid #ddd; border-radius:5px; margin-top:20px;'>
    <h3>Statistical Methods</h3>
    <ul>
        <li><strong>Global Structure:</strong> PERMANOVA (permutational MANOVA) with 999 permutations</li>
        <li><strong>Scoring Methods:</strong> GSVA (non-parametric) + Z-Score (parametric) with consensus reporting</li>
        <li><strong>Plasticity:</strong> Shannon entropy from signature scores
            <ul>
                <li>ANOVA (parametric): Assumes normality, more powerful when assumptions met</li>
                <li>Kruskal-Wallis (non-parametric): Robust to outliers, preferred for N<30 or non-normal data</li>
                <li>Shapiro-Wilk: Tests normality assumption (P<0.05 = non-normal)</li>
            </ul>
        </li>
        <li><strong>Trajectory Tests:</strong> 
            <ul>
                <li>Jonckheere-Terpstra: Tests for monotonic trend across ordered groups (more powerful than ANOVA for trajectories)</li>
                <li>Polynomial contrasts: Linear (steady progression) and Quadratic (spike/dip patterns)</li>
                <li>Consensus approach: Requires agreement between GSVA and Z-Score for robust findings</li>
            </ul>
        </li>
        <li><strong>Weighting Comparison:</strong> arrayWeights vs equal weights (automatic)</li>
        <li><strong>FDR Correction:</strong> Benjamini-Hochberg</li>
        <li><strong>Software:</strong> R ", as.character(R.version.string), ", limma ", as.character(packageVersion("limma")), 
        ", clinfun ", as.character(packageVersion("clinfun")), ", GSVA ", as.character(packageVersion("GSVA")), "</li>
    </ul>
    
    <h3>Interpretation Guide</h3>
    <ul>
        <li><strong>*** (P < 0.001):</strong> Highly significant</li>
        <li><strong>** (P < 0.01):</strong> Very significant</li>
        <li><strong>* (P < 0.05):</strong> Significant</li>
        <li><strong>. (P < 0.10):</strong> Trend</li>
        <li><strong>ns (P ≥ 0.10):</strong> Not significant</li>
    </ul>
    
    <h3>Key Concepts</h3>
    <ul>
        <li><strong>arrayWeights:</strong> Down-weights noisy samples. Can mask real biological switches in small datasets.</li>
        <li><strong>Jonckheere-Terpstra:</strong> Respects evolutionary order. Asks \"is there a monotonic trend?\" instead of \"are groups different?\"</li>
        <li><strong>Linear pattern:</strong> Signature steadily increases/decreases (e.g., EMT progression)</li>
        <li><strong>Quadratic pattern:</strong> Signature spikes or dips in middle stage (e.g., transient therapy response)</li>
        <li><strong>Consensus (ROBUST):</strong> Significant in BOTH GSVA AND Z-Score = highest confidence</li>
        <li><strong>Method-specific:</strong> Significant in only one method = moderate confidence, investigate further</li>
    </ul>
</div>

<p style='text-align:center; margin-top:20px; color:#666; font-size:12px;'>
Generated by <strong>v16.2 Omnibus Edition</strong> | ", as.character(Sys.time()), "
</p>
")

sink()
writeLines(capture.output(sessionInfo()), paste0(dirname(opt$out), "/sessionInfo.txt"))

cat("\n╔═══════════════════════════════════════════════╗\n")
cat("║   ANALYSIS COMPLETE - v16.2 OMNIBUS EDITION   ║\n")
cat("╚═══════════════════════════════════════════════╝\n\n")
cat("✓ Ran BOTH weighted and unweighted analyses\n")
cat("✓ Dual scoring: GSVA + Z-Score with consensus\n")
cat("✓ Fixed trajectory plot rendering (all lines show)\n")
cat("✓ Enhanced LLM summary with PCA, plasticity, correlations\n")
cat("✓ Jonckheere-Terpstra trajectory testing\n")
cat("✓ Polynomial contrasts (linear/quadratic)\n\n")
cat("📊 Key Findings:\n")
cat(sprintf("  • Consensus trends: %d/%d signatures\n", 
           sum(trajectory_results$Consensus_Trend == "ROBUST (Both)"), nrow(trajectory_results)))
cat(sprintf("  • Z-Score only: %d\n", sum(trajectory_results$Consensus_Trend == "Z-Score Only")))
cat(sprintf("  • GSVA only: %d\n", sum(trajectory_results$Consensus_Trend == "GSVA Only")))
cat(sprintf("  • Linear patterns: %d\n", sum(trajectory_results$Pattern == "Linear (monotonic)")))
cat(sprintf("  • Quadratic patterns: %d\n", sum(trajectory_results$Pattern == "Quadratic (spike/dip)")))
cat(sprintf("\n📝 Next: Review %s for biological interpretation!\n\n", basename(llm_file)))
