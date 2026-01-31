#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# GLOBAL SUBTYPES ANALYSIS v16.4 (ULTIMATE EDITION: v16.2 + v16.3 INTEGRATED)
# ------------------------------------------------------------------------------
# COMBINES:
# - v16.2: Enhanced LLM summary with PCA, plasticity, correlations, trajectory
# - v16.3: Comprehensive 2×4×N matrix (Scoring × Tests × Comparisons)
# - Professional significance heatmap
# - All correlation values reported
# - Complete pairwise testing framework
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
    library(clinfun)
    library(grid)
    library(gridExtra)
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
cat("LOG [1/13]: Loading Data...\n")

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
cat("LOG [2/13]: Global Structure Analysis...\n")

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
cat("LOG [3/13]: Scoring Signatures (DUAL: GSVA + Z-Score)...\n")

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

cat("  → Computing GSVA scores...\n")
gsva_res <- suppressWarnings(gsva(mat_sym, sigs, method="gsva", kcdf="Gaussian", verbose=FALSE))

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
# 4. COMPREHENSIVE PLASTICITY ANALYSIS (2×4 MATRIX)
# ==============================================================================
cat("LOG [4/13]: Comprehensive Plasticity Analysis (2×4 statistical matrix)...\n")

calc_entropy <- function(s) {
    s <- pmin(pmax(s, -10), 10)
    p <- exp(s)/sum(exp(s))
    -sum(p*log(p + 1e-10), na.rm=TRUE)
}

# Calculate plasticity for BOTH scoring methods
meta$Plasticity_ZScore <- apply(z_res, 2, calc_entropy)
meta$Plasticity_GSVA <- apply(gsva_res, 2, calc_entropy)

# Initialize results matrix
plasticity_matrix <- data.frame(
    Method = character(),
    Scoring = character(),
    Test = character(),
    Statistic = character(),
    P_value = numeric(),
    Significance = character(),
    stringsAsFactors = FALSE
)

# Function to run all tests for a given plasticity metric
run_plasticity_tests <- function(plasticity_col, scoring_name) {
    cat(sprintf("  → Testing %s plasticity...\n", scoring_name))

    # Normality test
    shapiro_test <- shapiro.test(plasticity_col)

    # ANOVA
    plast_aov <- aov(plasticity_col ~ meta$Classification)
    plast_summary <- summary(plast_aov)[[1]]
    anova_p <- plast_summary[["Pr(>F)"]][1]
    anova_f <- plast_summary[["F value"]][1]
    anova_df1 <- plast_summary[["Df"]][1]
    anova_df2 <- plast_summary[["Df"]][2]

    # Kruskal-Wallis
    kw_test <- kruskal.test(plasticity_col ~ meta$Classification)

    # Limma with weights
    design <- model.matrix(~0 + meta$Classification)
    colnames(design) <- levels(meta$Classification)
    mat_plast <- matrix(plasticity_col, nrow=1)
    colnames(mat_plast) <- names(plasticity_col)

    aw <- arrayWeights(mat_plast, design)
    fit_w <- lmFit(mat_plast, design, weights=aw)
    fit_w <- eBayes(fit_w)

    # Limma without weights
    fit_u <- lmFit(mat_plast, design)
    fit_u <- eBayes(fit_u)

    # Omnibus F-test from limma
    limma_w_p <- topTable(fit_w, number=1)$P.Value[1]
    limma_u_p <- topTable(fit_u, number=1)$P.Value[1]

    # Store results
    results <- rbind(
        data.frame(
            Method = "Plasticity",
            Scoring = scoring_name,
            Test = "Shapiro-Wilk (Normality)",
            Statistic = sprintf("W = %.3f", shapiro_test$statistic),
            P_value = shapiro_test$p.value,
            Significance = interpret_p(shapiro_test$p.value),
            stringsAsFactors = FALSE
        ),
        data.frame(
            Method = "Plasticity",
            Scoring = scoring_name,
            Test = "ANOVA",
            Statistic = sprintf("F(%d,%d) = %.3f", anova_df1, anova_df2, anova_f),
            P_value = anova_p,
            Significance = interpret_p(anova_p),
            stringsAsFactors = FALSE
        ),
        data.frame(
            Method = "Plasticity",
            Scoring = scoring_name,
            Test = "Kruskal-Wallis",
            Statistic = sprintf("H(%d) = %.3f", kw_test$parameter, kw_test$statistic),
            P_value = kw_test$p.value,
            Significance = interpret_p(kw_test$p.value),
            stringsAsFactors = FALSE
        ),
        data.frame(
            Method = "Plasticity",
            Scoring = scoring_name,
            Test = "Limma (Weighted)",
            Statistic = "Omnibus F",
            P_value = limma_w_p,
            Significance = interpret_p(limma_w_p),
            stringsAsFactors = FALSE
        ),
        data.frame(
            Method = "Plasticity",
            Scoring = scoring_name,
            Test = "Limma (Unweighted)",
            Statistic = "Omnibus F",
            P_value = limma_u_p,
            Significance = interpret_p(limma_u_p),
            stringsAsFactors = FALSE
        )
    )

    return(results)
}

# Run tests for both scoring methods
plasticity_matrix <- rbind(
    run_plasticity_tests(meta$Plasticity_ZScore, "Z-Score"),
    run_plasticity_tests(meta$Plasticity_GSVA, "GSVA")
)

# Export plasticity matrix
write.csv(plasticity_matrix, paste0(opt$out, "_Plasticity_Comprehensive_Matrix.csv"), row.names=FALSE)

# Add to stat_log
add_stat_log("Plasticity Tests Comprehensive", sprintf(
    "Ran 2×4 matrix of tests (2 scoring methods × 4 statistical tests)\n  %s",
    paste(capture.output(print(plasticity_matrix, row.names=FALSE)), collapse="\n  ")
))

# Plot plasticity for primary scoring method
meta$Plasticity <- meta$Plasticity_ZScore  # Use Z-Score as primary

# Get primary test recommendation
n_total <- nrow(meta)
shapiro_p_zscore <- plasticity_matrix$P_value[plasticity_matrix$Scoring == "Z-Score" & plasticity_matrix$Test == "Shapiro-Wilk (Normality)"]
if(n_total < 30 || shapiro_p_zscore < 0.05) {
    primary_test <- "Kruskal-Wallis"
    primary_p <- plasticity_matrix$P_value[plasticity_matrix$Scoring == "Z-Score" & plasticity_matrix$Test == "Kruskal-Wallis"]
} else {
    primary_test <- "ANOVA"
    primary_p <- plasticity_matrix$P_value[plasticity_matrix$Scoring == "Z-Score" & plasticity_matrix$Test == "ANOVA"]
}

p_plast <- ggplot(meta, aes(x=Classification, y=Plasticity, fill=Classification)) +
    geom_violin(alpha=0.3, trim=FALSE) +
    geom_boxplot(width=0.2, fill="white", outlier.shape=NA) +
    geom_jitter(width=0.1, size=3, alpha=0.7) +
    stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="white") +
    scale_fill_manual(values=GROUP_COLORS) +
    labs(title="Cellular Plasticity (Shannon Entropy - Z-Score)",
         subtitle=sprintf("Recommended: %s P=%s %s | See matrix for all tests", 
                         primary_test, safe_format(primary_p), interpret_p(primary_p)),
         caption="Diamond = mean, box = median ± IQR, dots = individual samples") +
    theme_publication() +
    theme(legend.position="none")

ggsave(paste0(opt$out, "_Plasticity_mqc.png"), p_plast, width=7, height=6)

# Descriptive statistics
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

# ==============================================================================
# 5. DUAL-WEIGHTING DIFFERENTIAL ANALYSIS
# ==============================================================================
cat("LOG [5/13]: Dual-Weighting Differential Analysis...\n")

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

    levs <- levels(metadata$Classification)
    contrast_formulas <- c()
    contrast_names <- c()

    # Generate ALL pairwise contrasts
    for(i in 1:(length(levs)-1)) {
        for(j in (i+1):length(levs)) {
            contrast_formulas <- c(contrast_formulas, sprintf("%s - %s", levs[j], levs[i]))
            contrast_names <- c(contrast_names, sprintf("%s_vs_%s",
                                                        gsub("[^A-Za-z0-9]", "", levs[j]),
                                                        gsub("[^A-Za-z0-9]", "", levs[i])))
        }
    }

    cont.matrix <- makeContrasts(contrasts=contrast_formulas, levels=design)
    colnames(cont.matrix) <- contrast_names

    # Polynomial contrasts for trajectory
    design_poly <- model.matrix(~ poly(as.numeric(metadata$Classification), 2))
    colnames(design_poly) <- c("Intercept", "Linear", "Quadratic")

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

cat("\n╔═══════════════════════════════════════════════╗\n")
cat("║   ANALYSIS 1: WITH arrayWeights (Standard)   ║\n")
cat("╚═══════════════════════════════════════════════╝\n")
res_weighted <- run_limma_analysis(final_scores, meta, use_weights=TRUE, label="WEIGHTED")

cat("\n╔═══════════════════════════════════════════════╗\n")
cat("║  ANALYSIS 2: WITHOUT arrayWeights (Unbiased)  ║\n")
cat("╚═══════════════════════════════════════════════╝\n")
res_unweighted <- run_limma_analysis(final_scores, meta, use_weights=FALSE, label="UNWEIGHTED")

all_coefs <- res_weighted$contrasts
levs <- res_weighted$levs

# ==============================================================================
# 6. COMPREHENSIVE TRAJECTORY TESTING (2×4×N MATRIX)
# ==============================================================================
cat("\nLOG [6/13]: Comprehensive Trajectory Testing (JT + ALL Pairwise)...\n")

stage_numeric <- as.numeric(meta$Classification)

# Initialize comprehensive results dataframe
trajectory_results <- data.frame(
    Signature = rownames(final_scores),
    stringsAsFactors = FALSE
)

# Generate all pairwise combinations
pairwise_combos <- list()
for(i in 1:(length(levs)-1)) {
    for(j in (i+1):length(levs)) {
        combo_name <- paste0(gsub("[^A-Za-z0-9]", "", levs[j]), "_vs_", gsub("[^A-Za-z0-9]", "", levs[i]))
        pairwise_combos[[combo_name]] <- c(levs[i], levs[j])
    }
}

cat(sprintf("  → Testing %d pairwise comparisons: %s\n",
            length(pairwise_combos), paste(names(pairwise_combos), collapse=", ")))

# Function to run comprehensive tests for one scoring method
run_comprehensive_trajectory <- function(curr_scores, method_name) {
    cat(sprintf("\n  → Running comprehensive tests for %s...\n", method_name))

    results_list <- list()

    # 1. Global Jonckheere-Terpstra (monotonic trend across all groups)
    cat("    • Global JT test (all groups)...\n")
    for(i in 1:nrow(curr_scores)) {
        sig_name <- rownames(curr_scores)[i]
        scores <- curr_scores[i, ]

        jt <- tryCatch({
            jonckheere.test(scores, stage_numeric, nperm=1000)
        }, error = function(e) {
            list(statistic = NA, p.value = NA)
        })

        n <- length(scores)
        g_sizes <- table(stage_numeric)
        expected_mean <- (n^2 - sum(g_sizes^2)) / 4

        results_list[[sig_name]][[paste0(method_name, "_Global_JT_Stat")]] <- jt$statistic
        results_list[[sig_name]][[paste0(method_name, "_Global_JT_P")]] <- jt$p.value
        results_list[[sig_name]][[paste0(method_name, "_Global_JT_Direction")]] <- ifelse(
            is.na(jt$statistic), "Unknown",
            ifelse(jt$statistic > expected_mean, "Increasing", "Decreasing")
        )
    }

    # 2. Pairwise tests for EACH combination
    for(combo_name in names(pairwise_combos)) {
        groups_to_test <- pairwise_combos[[combo_name]]
        cat(sprintf("    • Pairwise: %s...\n", combo_name))

        # Subset data
        subset_idx <- meta$Classification %in% groups_to_test
        subset_scores <- curr_scores[, subset_idx, drop=FALSE]
        subset_meta <- meta[subset_idx, , drop=FALSE]
        subset_meta$Classification <- droplevels(subset_meta$Classification)

        for(i in 1:nrow(subset_scores)) {
            sig_name <- rownames(subset_scores)[i]

            # T-test
            group1_vals <- subset_scores[i, subset_meta$Classification == groups_to_test[1]]
            group2_vals <- subset_scores[i, subset_meta$Classification == groups_to_test[2]]

            ttest <- tryCatch({
                t.test(group2_vals, group1_vals)
            }, error = function(e) {
                list(statistic = NA, p.value = NA, estimate = c(NA, NA))
            })

            # Wilcoxon test
            wilcox <- tryCatch({
                wilcox.test(group2_vals, group1_vals)
            }, error = function(e) {
                list(statistic = NA, p.value = NA)
            })

            # Limma (weighted)
            design_pw <- model.matrix(~0 + subset_meta$Classification)
            colnames(design_pw) <- levels(subset_meta$Classification)
            cont_pw <- makeContrasts(
                contrasts = sprintf("%s - %s", groups_to_test[2], groups_to_test[1]),
                levels = design_pw
            )

            aw_pw <- arrayWeights(matrix(subset_scores[i,], nrow=1), design_pw)
            fit_pw_w <- lmFit(matrix(subset_scores[i,], nrow=1), design_pw, weights=aw_pw)
            fit_pw_w <- contrasts.fit(fit_pw_w, cont_pw)
            fit_pw_w <- eBayes(fit_pw_w)
            tt_pw_w <- topTable(fit_pw_w, number=1)

            # Limma (unweighted)
            fit_pw_u <- lmFit(matrix(subset_scores[i,], nrow=1), design_pw)
            fit_pw_u <- contrasts.fit(fit_pw_u, cont_pw)
            fit_pw_u <- eBayes(fit_pw_u)
            tt_pw_u <- topTable(fit_pw_u, number=1)

            # Store results
            results_list[[sig_name]][[paste0(method_name, "_", combo_name, "_TTest_P")]] <- ttest$p.value
            results_list[[sig_name]][[paste0(method_name, "_", combo_name, "_TTest_T")]] <- ttest$statistic
            results_list[[sig_name]][[paste0(method_name, "_", combo_name, "_Wilcox_P")]] <- wilcox$p.value
            results_list[[sig_name]][[paste0(method_name, "_", combo_name, "_Limma_W_P")]] <- tt_pw_w$adj.P.Val[1]
            results_list[[sig_name]][[paste0(method_name, "_", combo_name, "_Limma_W_logFC")]] <- tt_pw_w$logFC[1]
            results_list[[sig_name]][[paste0(method_name, "_", combo_name, "_Limma_U_P")]] <- tt_pw_u$adj.P.Val[1]
            results_list[[sig_name]][[paste0(method_name, "_", combo_name, "_Limma_U_logFC")]] <- tt_pw_u$logFC[1]
        }
    }

    # Convert results list to dataframe
    results_df <- do.call(rbind, lapply(names(results_list), function(sig) {
        data.frame(Signature = sig, results_list[[sig]], stringsAsFactors = FALSE)
    }))

    return(results_df)
}

# Run for BOTH scoring methods
zscore_traj <- run_comprehensive_trajectory(z_res, "ZScore")
gsva_traj <- run_comprehensive_trajectory(gsva_res, "GSVA")

# Merge results
trajectory_results <- merge(trajectory_results, zscore_traj, by="Signature", all=TRUE)
trajectory_results <- merge(trajectory_results, gsva_traj, by="Signature", all=TRUE)

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
    jt_p <- trajectory_results$ZScore_Global_JT_P[i]
    
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
    z_sig <- !is.na(trajectory_results$ZScore_Global_JT_P[i]) && trajectory_results$ZScore_Global_JT_P[i] < 0.05
    g_sig <- !is.na(trajectory_results$GSVA_Global_JT_P[i]) && trajectory_results$GSVA_Global_JT_P[i] < 0.05
    
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

# Export comprehensive trajectory results
write.csv(trajectory_results, paste0(opt$out, "_Trajectory_Comprehensive_ALL_Tests.csv"), row.names=FALSE)

cat(sprintf("\n✓ Comprehensive trajectory testing complete: %d signatures × %d comparisons × 2 methods × 4 tests\n",
           nrow(trajectory_results), length(pairwise_combos) + 1))

# ==============================================================================
# 7. CREATE SIGNIFICANCE MATRIX HEATMAP
# ==============================================================================
cat("\nLOG [7/13]: Creating Significance Matrix Heatmap...\n")

# Build comprehensive significance matrix
create_sig_matrix <- function() {

    sig_matrix <- data.frame(
        Signature = character(),
        Comparison = character(),
        Scoring = character(),
        Test = character(),
        P_value = numeric(),
        Is_Sig = character(),
        stringsAsFactors = FALSE
    )

    # Process each signature
    for(sig in rownames(final_scores)) {

        # Global JT tests
        for(scoring in c("ZScore", "GSVA")) {
            jt_col <- paste0(scoring, "_Global_JT_P")
            if(jt_col %in% colnames(trajectory_results)) {
                p_val <- trajectory_results[trajectory_results$Signature == sig, jt_col]
                sig_matrix <- rbind(sig_matrix, data.frame(
                    Signature = sig,
                    Comparison = "Global_Trajectory",
                    Scoring = scoring,
                    Test = "JT",
                    P_value = p_val,
                    Is_Sig = ifelse(is.na(p_val), "NA", ifelse(p_val < 0.05, "✓", "✗")),
                    stringsAsFactors = FALSE
                ))
            }
        }

        # Pairwise tests
        for(combo_name in names(pairwise_combos)) {
            for(scoring in c("ZScore", "GSVA")) {

                # T-test
                ttest_col <- paste0(scoring, "_", combo_name, "_TTest_P")
                if(ttest_col %in% colnames(trajectory_results)) {
                    p_val <- trajectory_results[trajectory_results$Signature == sig, ttest_col]
                    sig_matrix <- rbind(sig_matrix, data.frame(
                        Signature = sig,
                        Comparison = combo_name,
                        Scoring = scoring,
                        Test = "T-test",
                        P_value = p_val,
                        Is_Sig = ifelse(is.na(p_val), "NA", ifelse(p_val < 0.05, "✓", "✗")),
                        stringsAsFactors = FALSE
                    ))
                }

                # Wilcoxon
                wilcox_col <- paste0(scoring, "_", combo_name, "_Wilcox_P")
                if(wilcox_col %in% colnames(trajectory_results)) {
                    p_val <- trajectory_results[trajectory_results$Signature == sig, wilcox_col]
                    sig_matrix <- rbind(sig_matrix, data.frame(
                        Signature = sig,
                        Comparison = combo_name,
                        Scoring = scoring,
                        Test = "Wilcoxon",
                        P_value = p_val,
                        Is_Sig = ifelse(is.na(p_val), "NA", ifelse(p_val < 0.05, "✓", "✗")),
                        stringsAsFactors = FALSE
                    ))
                }

                # Limma Weighted
                limma_w_col <- paste0(scoring, "_", combo_name, "_Limma_W_P")
                if(limma_w_col %in% colnames(trajectory_results)) {
                    p_val <- trajectory_results[trajectory_results$Signature == sig, limma_w_col]
                    sig_matrix <- rbind(sig_matrix, data.frame(
                        Signature = sig,
                        Comparison = combo_name,
                        Scoring = scoring,
                        Test = "Limma_W",
                        P_value = p_val,
                        Is_Sig = ifelse(is.na(p_val), "NA", ifelse(p_val < 0.05, "✓", "✗")),
                        stringsAsFactors = FALSE
                    ))
                }

                # Limma Unweighted
                limma_u_col <- paste0(scoring, "_", combo_name, "_Limma_U_P")
                if(limma_u_col %in% colnames(trajectory_results)) {
                    p_val <- trajectory_results[trajectory_results$Signature == sig, limma_u_col]
                    sig_matrix <- rbind(sig_matrix, data.frame(
                        Signature = sig,
                        Comparison = combo_name,
                        Scoring = scoring,
                        Test = "Limma_U",
                        P_value = p_val,
                        Is_Sig = ifelse(is.na(p_val), "NA", ifelse(p_val < 0.05, "✓", "✗")),
                        stringsAsFactors = FALSE
                    ))
                }
            }
        }
    }

    return(sig_matrix)
}

sig_matrix <- create_sig_matrix()
write.csv(sig_matrix, paste0(opt$out, "_Significance_Matrix_Full.csv"), row.names=FALSE)

# Create visual heatmap
sig_matrix$Test_Label <- paste0(sig_matrix$Scoring, "_", sig_matrix$Test)
sig_matrix$Comp_Test <- paste0(sig_matrix$Comparison, "\n", sig_matrix$Test_Label)

# Prepare matrix for heatmap
plot_matrix <- sig_matrix %>%
    dplyr::select(Signature, Comp_Test, Is_Sig) %>%
    pivot_wider(names_from = Comp_Test, values_from = Is_Sig, values_fill = "NA")

plot_matrix_mat <- as.matrix(plot_matrix[, -1])
rownames(plot_matrix_mat) <- plot_matrix$Signature

# Create color mapping
color_map <- c("✓" = "#2ecc71", "✗" = "#e74c3c", "NA" = "#95a5a6")

png(paste0(opt$out, "_Significance_Matrix_Heatmap_mqc.png"), width=16, height=10, units="in", res=300)

par(mar=c(12, 15, 4, 2))
n_rows <- nrow(plot_matrix_mat)
n_cols <- ncol(plot_matrix_mat)

plot(1, type="n", xlim=c(0.5, n_cols+0.5), ylim=c(0.5, n_rows+0.5),
     xlab="", ylab="", xaxt="n", yaxt="n", main="Comprehensive Significance Matrix (2×4×N Testing)")

# Draw cells
for(i in 1:n_rows) {
    for(j in 1:n_cols) {
        val <- plot_matrix_mat[i, j]
        rect(j-0.5, i-0.5, j+0.5, i+0.5, col=color_map[val], border="white", lwd=2)
        text(j, i, val, cex=1.2, font=2)
    }
}

axis(1, at=1:n_cols, labels=colnames(plot_matrix_mat), las=2, cex.axis=0.7)
axis(2, at=1:n_rows, labels=rownames(plot_matrix_mat), las=2, cex.axis=0.9)

legend("topright", legend=c("Significant (P<0.05)", "Not Significant", "Not Available"),
       fill=c("#2ecc71", "#e74c3c", "#95a5a6"), border="black", cex=0.9)

dev.off()

cat("✓ Significance matrix heatmap created\n")

# ==============================================================================
# 8. UNIFIED TRAJECTORY PLOT
# ==============================================================================
cat("\nLOG [8/13]: Trajectory Visualization...\n")

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

traj_summary <- traj_data %>%
    group_by(Signature, Class, Stage) %>%
    summarise(Mean = mean(Score, na.rm=TRUE), SE = sd(Score, na.rm=TRUE)/sqrt(n()), .groups="drop")

# Order by JT statistic
if("ZScore_Global_JT_Stat" %in% colnames(trajectory_results)) {
    trend_order <- trajectory_results %>%
        arrange(desc(ZScore_Global_JT_Stat)) %>%
        pull(Signature)
} else {
    trend_order <- rownames(final_scores)
}

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
         subtitle="Ordered by JT trend statistic | Lines = group means ± SE",
         x="Stage", y="Z-Score Expression") +
    theme_publication(base_size=11) +
    theme(legend.position = "bottom")

ggsave(paste0(opt$out, "_Unified_Trajectories_mqc.png"), p_traj, width=14, height=12)

# ==============================================================================
# 9. CORRELATION ANALYSIS
# ==============================================================================
cat("LOG [9/13]: Correlation Analysis...\n")

sig_cor <- cor(t(final_scores), method="pearson")

# Export FULL correlation matrix
write.csv(sig_cor, paste0(opt$out, "_Signature_Correlation_Matrix_FULL.csv"))

cor_pairs <- which(abs(sig_cor) > 0.6 & upper.tri(sig_cor), arr.ind=TRUE)
if(nrow(cor_pairs) > 0) {
    cor_summary <- data.frame(
        Sig1 = rownames(sig_cor)[cor_pairs[,1]],
        Sig2 = rownames(sig_cor)[cor_pairs[,2]],
        Correlation = sig_cor[cor_pairs]
    )
    cor_summary <- cor_summary[order(abs(cor_summary$Correlation), decreasing=TRUE), ]
    write.csv(cor_summary, paste0(opt$out, "_Signature_Correlations_High.csv"), row.names=FALSE)
    
    add_stat_log("Signature Co-evolution (High Correlations)", sprintf(
        "Detected %d high correlations (|r| > 0.6):\n  %s",
        nrow(cor_summary),
        paste(capture.output(print(head(cor_summary, 10), row.names=FALSE)), collapse="\n  ")
    ))
} else {
    cor_summary <- data.frame()
    add_stat_log("Signature Co-evolution", "No strong correlations (|r| > 0.6) detected")
}

# ==============================================================================
# 10. HEATMAPS
# ==============================================================================
cat("LOG [10/13]: Generating Heatmaps...\n")

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
# 11. ENHANCED LLM SUMMARY (ULTIMATE EDITION)
# ==============================================================================
cat("LOG [11/13]: Generating Ultimate LLM Summary...\n")

llm_summary <- paste0(
    "═══════════════════════════════════════════════════════════════════\n",
    "GBM LITT THERAPY SUBTYPE EVOLUTION - ULTIMATE ANALYSIS v16.4\n",
    "═══════════════════════════════════════════════════════════════════\n\n",
    "EXPERIMENTAL DESIGN:\n",
    sprintf("  Trajectory: %s\n", paste(levs, collapse=" → ")),
    sprintf("  Samples: %s\n", paste(names(group_sizes), "=", group_sizes, collapse=", ")),
    sprintf("  Total signatures tested: %d\n", length(sigs)),
    sprintf("  Statistical framework: 2×4×%d matrix (Scoring×Tests×Comparisons)\n", length(pairwise_combos)+1),
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
    "2. CELLULAR PLASTICITY (Comprehensive 2×4 Matrix)\n",
    "═══════════════════════════════════════════════════════════════════\n"
)

# Add plasticity results
for(scoring in c("Z-Score", "GSVA")) {
    llm_summary <- paste0(llm_summary, sprintf("\n%s Plasticity:\n", scoring))
    scoring_results <- plasticity_matrix[plasticity_matrix$Scoring == scoring, ]
    for(i in 1:nrow(scoring_results)) {
        llm_summary <- paste0(llm_summary,
            sprintf("  %s: %s, P = %s %s\n",
                   scoring_results$Test[i],
                   scoring_results$Statistic[i],
                   safe_format(scoring_results$P_value[i]),
                   scoring_results$Significance[i]))
    }
}

llm_summary <- paste0(llm_summary, "\n  Group Statistics (Z-Score):\n")
for(i in 1:nrow(plast_stats)) {
    llm_summary <- paste0(llm_summary,
        sprintf("    %s: Mean=%.3f±%.3f, Median=%.3f, IQR=%.3f (n=%d)\n",
                plast_stats$Classification[i], plast_stats$Mean[i], plast_stats$SD[i],
                plast_stats$Median[i], plast_stats$IQR[i], plast_stats$N[i]))
}

llm_summary <- paste0(llm_summary, "\n",
    "═══════════════════════════════════════════════════════════════════\n",
    "3. SIGNATURE CO-EVOLUTION (Complete Correlation Matrix)\n",
    "═══════════════════════════════════════════════════════════════════\n"
)

# Add ALL pairwise correlations
llm_summary <- paste0(llm_summary, sprintf("  Total signature pairs analyzed: %d\n\n", (nrow(sig_cor)*(nrow(sig_cor)-1))/2))

if(nrow(cor_summary) > 0) {
    llm_summary <- paste0(llm_summary, sprintf("  High Correlations (|r| > 0.6): %d pairs\n\n", nrow(cor_summary)))
    for(i in 1:min(15, nrow(cor_summary))) {
        llm_summary <- paste0(llm_summary,
            sprintf("    %s ↔ %s: r=%.3f\n",
                    cor_summary$Sig1[i], cor_summary$Sig2[i], cor_summary$Correlation[i]))
    }
    if(nrow(cor_summary) > 15) {
        llm_summary <- paste0(llm_summary, sprintf("\n    ... and %d more (see full matrix CSV)\n", nrow(cor_summary) - 15))
    }
} else {
    llm_summary <- paste0(llm_summary, "  No strong correlations (|r| > 0.6) detected\n")
}

llm_summary <- paste0(llm_summary, "\n",
    "═══════════════════════════════════════════════════════════════════\n",
    "4. TRAJECTORY TRENDS (Dual Method Consensus + Patterns)\n",
    "═══════════════════════════════════════════════════════════════════\n"
)

# Add trajectory pattern summary
llm_summary <- paste0(llm_summary, "  Pattern Distribution:\n")
pattern_table <- table(trajectory_results$Pattern)
for(pattern in names(pattern_table)) {
    llm_summary <- paste0(llm_summary, sprintf("    %s: %d signatures\n", pattern, pattern_table[pattern]))
}

llm_summary <- paste0(llm_summary, "\n  Consensus Trends:\n")
consensus_table <- table(trajectory_results$Consensus_Trend)
for(cons in names(consensus_table)) {
    llm_summary <- paste0(llm_summary, sprintf("    %s: %d signatures\n", cons, consensus_table[cons]))
}

sig_trends <- trajectory_results[trajectory_results$Consensus_Trend != "None", ]
if(nrow(sig_trends) > 0) {
    llm_summary <- paste0(llm_summary, "\n  Significant Monotonic Trends:\n")
    sig_trends <- sig_trends[order(match(sig_trends$Consensus_Trend, c("ROBUST (Both)", "Z-Score Only", "GSVA Only"))), ]
    for(i in 1:nrow(sig_trends)) {
        sig <- sig_trends$Signature[i]
        llm_summary <- paste0(llm_summary, sprintf(
            "    %s [%s]:\n      Z-Score JT: P=%s (Stat=%.1f, Direction=%s)\n      GSVA JT: P=%s (Stat=%.1f, Direction=%s)\n      Pattern: %s\n",
            sig, sig_trends$Consensus_Trend[i],
            safe_format(sig_trends$ZScore_Global_JT_P[i]), sig_trends$ZScore_Global_JT_Stat[i], sig_trends$ZScore_Global_JT_Direction[i],
            safe_format(sig_trends$GSVA_Global_JT_P[i]), sig_trends$GSVA_Global_JT_Stat[i], sig_trends$GSVA_Global_JT_Direction[i],
            sig_trends$Pattern[i]
        ))
    }
}

llm_summary <- paste0(llm_summary, "\n",
    "═══════════════════════════════════════════════════════════════════\n",
    "5. PAIRWISE COMPARISONS (2×4 Matrix per Comparison)\n",
    "═══════════════════════════════════════════════════════════════════\n",
    sprintf("  Total comparisons tested: %d (Global + %d Pairwise)\n",
            length(pairwise_combos)+1, length(pairwise_combos)),
    sprintf("  Tests per comparison: 8 (2 scoring × 4 statistical tests)\n\n")
)

# Add pairwise results
for(combo_name in names(pairwise_combos)) {
    llm_summary <- paste0(llm_summary, sprintf("  %s (%s):\n", 
                                                combo_name, 
                                                paste(pairwise_combos[[combo_name]], collapse=" vs ")))

    # Get significant signatures for this comparison
    combo_sigs <- sig_matrix[sig_matrix$Comparison == combo_name & sig_matrix$Is_Sig == "✓", ]
    
    if(nrow(combo_sigs) > 0) {
        # Group by signature and count significant tests
        sig_counts <- combo_sigs %>%
            group_by(Signature) %>%
            summarise(N_Sig_Tests = n(), .groups="drop") %>%
            arrange(desc(N_Sig_Tests))
        
        top_sigs <- head(sig_counts, 5)
        
        for(j in 1:nrow(top_sigs)) {
            sig_name <- top_sigs$Signature[j]
            n_sig <- top_sigs$N_Sig_Tests[j]
            
            # Get all tests for this signature
            sig_tests <- combo_sigs[combo_sigs$Signature == sig_name, ]
            
            llm_summary <- paste0(llm_summary, sprintf("    • %s (%d/8 tests significant):\n", sig_name, n_sig))
            for(k in 1:nrow(sig_tests)) {
                llm_summary <- paste0(llm_summary,
                    sprintf("      %s (%s): P = %s\n",
                           sig_tests$Test[k],
                           sig_tests$Scoring[k],
                           safe_format(sig_tests$P_value[k])))
            }
        }
        
        if(nrow(sig_counts) > 5) {
            llm_summary <- paste0(llm_summary, sprintf("    ... and %d more signatures (see matrix)\n", nrow(sig_counts) - 5))
        }
    } else {
        llm_summary <- paste0(llm_summary, "    No significant changes detected in any test\n")
    }
    llm_summary <- paste0(llm_summary, "\n")
}

llm_summary <- paste0(llm_summary,
    "═══════════════════════════════════════════════════════════════════\n",
    "INTERPRETATION GUIDE FOR LLM\n",
    "═══════════════════════════════════════════════════════════════════\n",
    "Statistical Framework:\n",
    "  Dimension 1 (Scoring): GSVA vs Z-Score\n",
    "  Dimension 2 (Tests): T-test, Wilcoxon, Limma-W, Limma-U\n",
    "  Dimension 3 (Comparisons): Global trajectory + All pairwise\n",
    sprintf("  Total tests per signature: %d\n\n", 2*4*(length(pairwise_combos)+1)),
    
    "Consensus Interpretation:\n",
    "  ROBUST (Both): Significant in BOTH Z-Score AND GSVA → Highest confidence\n",
    "  Method-specific: Significant in only one method → Moderate confidence\n",
    "  Agreement across tests: Multiple tests significant → Robust finding\n",
    "  Weighted vs Unweighted: Differences indicate sample quality effects\n\n",
    
    "Pattern Classification:\n",
    "  Linear (monotonic): Steady progression (e.g., EMT evolution)\n",
    "  Quadratic (spike/dip): Transient state (e.g., therapy shock)\n",
    "  Weak trend: Suggestive but not definitive\n",
    "  No trend: Stable across trajectory\n\n",
    
    "Clinical Relevance Priority:\n",
    "  1. Primary vs Recurrent = Core therapy effect\n",
    "  2. Culture vs Primary = Microenvironment influence\n",
    "  3. Culture vs Recurrent = Total evolution span\n",
    "  4. Global trajectory = Overall pattern\n\n",
    
    "Focus on signatures with:\n",
    "  - ROBUST consensus (both scoring methods)\n",
    "  - Multiple significant tests (≥3/4)\n",
    "  - Biological relevance to GBM recurrence\n",
    "  - Clear directional trend (increasing/decreasing)\n\n",
    
    "Significance Codes:\n",
    "  *** P < 0.001 (highly significant)\n",
    "  **  P < 0.01  (very significant)\n",
    "  *   P < 0.05  (significant)\n",
    "  .   P < 0.10  (trend)\n",
    "  ns  P ≥ 0.10  (not significant)\n",
    "═══════════════════════════════════════════════════════════════════\n"
)

llm_file <- paste0(opt$out, "_llm_summary_ULTIMATE.txt")
writeLines(llm_summary, llm_file)
cat(sprintf("\n✓ Ultimate LLM summary: %s\n", basename(llm_file)))

# ==============================================================================
# 12. HTML REPORT
# ==============================================================================
cat("LOG [12/13]: Generating HTML Report...\n")

summary_html <- paste0(dirname(opt$out), "/analysis_summary_mqc.html")
sink(summary_html)

cat("
<div style='background-color:#f8f9fa; padding:20px; border:2px solid #667eea; border-radius:8px;'>
    <h2 style='color:#667eea;'>🧬 Litt Therapy Subtype Evolution - ULTIMATE ANALYSIS v16.4</h2>

    <div style='background:#d1ecf1; border-left:4px solid#0c5460; padding:10px; margin:15px 0;'>
    <strong>🆕 v16.4 ULTIMATE FEATURES:</strong><br>
    ✓ Complete integration of v16.2 enhanced LLM reporting<br>
    ✓ Complete integration of v16.3 comprehensive 2×4×N matrix<br>
    ✓ Professional significance heatmap with checkmarks (✓/✗)<br>
    ✓ All correlation values in LLM summary<br>
    ✓ PCA gene drivers, plasticity statistics, trajectory patterns<br>
    ✓ Full pairwise testing framework
    </div>

    <h3>Analysis Framework</h3>
    <table style='border-collapse:collapse; width:100%; margin:15px 0;'>
        <tr style='background:#667eea; color:white;'>
            <th style='padding:8px; border:1px solid #555;'>Component</th>
            <th style='padding:8px; border:1px solid #555;'>Details</th>
        </tr>
        <tr><td style='padding:8px; border:1px solid #ddd;'><strong>Scoring Methods</strong></td>
            <td style='padding:8px; border:1px solid #ddd;'>GSVA + Z-Score (dual approach)</td></tr>
        <tr><td style='padding:8px; border:1px solid #ddd;'><strong>Plasticity Tests</strong></td>
            <td style='padding:8px; border:1px solid #ddd;'>Shapiro-Wilk, ANOVA, Kruskal-Wallis, Limma-W, Limma-U</td></tr>
        <tr><td style='padding:8px; border:1px solid #ddd;'><strong>Trajectory Tests</strong></td>
            <td style='padding:8px; border:1px solid #ddd;'>Global JT + ", length(pairwise_combos), " Pairwise × 4 tests</td></tr>
        <tr><td style='padding:8px; border:1px solid #ddd;'><strong>Total Tests/Signature</strong></td>
            <td style='padding:8px; border:1px solid #ddd;'><strong>", 2*4*(length(pairwise_combos)+1), "</strong></td></tr>
    </table>

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

    <h3>Plasticity Analysis (2×4 Matrix)</h3>
    <p>All ", nrow(plasticity_matrix), " tests shown in comprehensive matrix CSV</p>
    <p><strong>Recommendation:</strong> ", 
    if(n_total < 30 || shapiro_p_zscore < 0.05) {
        sprintf("Use <strong>Kruskal-Wallis</strong> (N=%d, data is %s)", 
                n_total, if(shapiro_p_zscore < 0.05) "non-normal" else "small sample")
    } else {
        sprintf("Use <strong>ANOVA</strong> (N=%d, data is normally distributed)", n_total)
    },
    "</p>

    <h3>Trajectory Pattern Summary</h3>
    <table style='border-collapse:collapse; width:100%;'>
        <tr style='background:#667eea; color:white;'>
            <th style='padding:8px;'>Pattern</th>
            <th style='padding:8px;'>Count</th>
        </tr>
")

pattern_table <- table(trajectory_results$Pattern)
for(pattern in names(pattern_table)) {
    cat(sprintf("<tr><td style='padding:8px;'>%s</td><td style='padding:8px;'>%d</td></tr>\n",
               pattern, pattern_table[pattern]))
}

cat("    </table>

    <h3>Consensus Trends</h3>
    <table style='border-collapse:collapse; width:100%;'>
        <tr style='background:#667eea; color:white;'>
            <th style='padding:8px;'>Consensus</th>
            <th style='padding:8px;'>Count</th>
        </tr>
")

consensus_table <- table(trajectory_results$Consensus_Trend)
for(cons in names(consensus_table)) {
    cat(sprintf("<tr><td style='padding:8px;'>%s</td><td style='padding:8px;'>%d</td></tr>\n",
               cons, consensus_table[cons]))
}

cat("    </table>

    <h3>Key Pairwise Comparisons</h3>
    <table style='border-collapse:collapse; width:100%;'>
        <tr style='background:#667eea; color:white;'>
            <th style='padding:8px;'>Comparison</th>
            <th style='padding:8px;'>Groups</th>
        </tr>
")

for(combo_name in names(pairwise_combos)) {
    cat(sprintf("        <tr><td style='padding:8px; border:1px solid #ddd;'>%s</td>
            <td style='padding:8px; border:1px solid #ddd;'>%s</td></tr>\n",
            combo_name, paste(pairwise_combos[[combo_name]], collapse=" vs ")))
}

cat("    </table>
</div>

<div style='background:#e7f3ff; padding:20px; border:2px solid #667eea; border-radius:8px; margin-top:20px;'>
    <h2 style='color:#667eea;'>📄 Ultimate LLM Interpretation Summary</h2>
    <p style='color:#666; font-size:13px;'>
    <strong>Complete Analysis Summary:</strong> Includes PCA components, plasticity statistics, 
    all correlations, trajectory patterns, and comprehensive pairwise testing results.<br>
    <strong>File:</strong> <code>", basename(llm_file), "</code><br>
    <strong>Purpose:</strong> Copy into ChatGPT/Claude/Gemini for biological interpretation
    </p>
</div>

<div style='background:#fff; padding:15px; border:1px solid #ddd; border-radius:5px; margin-top:20px;'>
    <h3>📊 Generated Files</h3>
    <ul>
        <li><code>", basename(paste0(opt$out, "_Significance_Matrix_Heatmap_mqc.png")), "</code> - Visual significance matrix</li>
        <li><code>", basename(paste0(opt$out, "_Significance_Matrix_Full.csv")), "</code> - Complete test results</li>
        <li><code>", basename(paste0(opt$out, "_Trajectory_Comprehensive_ALL_Tests.csv")), "</code> - All trajectory tests</li>
        <li><code>", basename(paste0(opt$out, "_Plasticity_Comprehensive_Matrix.csv")), "</code> - Plasticity 2×4 matrix</li>
        <li><code>", basename(paste0(opt$out, "_Signature_Correlation_Matrix_FULL.csv")), "</code> - Complete correlation matrix</li>
        <li><code>", basename(llm_file), "</code> - Ultimate LLM summary</li>
    </ul>
</div>

<div style='background:#fff; padding:15px; border:1px solid #ddd; border-radius:5px; margin-top:20px;'>
    <h3>Statistical Methods Summary</h3>
    <ul>
        <li><strong>Global Structure:</strong> PERMANOVA (999 permutations), PCA with gene loadings</li>
        <li><strong>Scoring:</strong> GSVA (non-parametric) + Z-Score (parametric) with correlation assessment</li>
        <li><strong>Plasticity:</strong> Shannon entropy with 2×4 testing matrix
            <ul>
                <li>Normality: Shapiro-Wilk</li>
                <li>Parametric: ANOVA, Limma (weighted/unweighted)</li>
                <li>Non-parametric: Kruskal-Wallis</li>
            </ul>
        </li>
        <li><strong>Trajectory:</strong> 
            <ul>
                <li>Global monotonic trend: Jonckheere-Terpstra (1000 permutations)</li>
                <li>Pairwise comparisons: T-test, Wilcoxon, Limma-W, Limma-U</li>
                <li>Pattern classification: Linear vs Quadratic (polynomial contrasts)</li>
            </ul>
        </li>
        <li><strong>Multiple Testing:</strong> Benjamini-Hochberg FDR correction</li>
        <li><strong>Software:</strong> R ", as.character(R.version.string), 
        ", limma ", as.character(packageVersion("limma")), 
        ", clinfun ", as.character(packageVersion("clinfun")), 
        ", GSVA ", as.character(packageVersion("GSVA")), "</li>
    </ul>
</div>

<p style='text-align:center; margin-top:20px; color:#666; font-size:12px;'>
Generated by <strong>v16.4 Ultimate Edition</strong> | ", as.character(Sys.time()), "
</p>
")

sink()

# ==============================================================================
# 13. EXPORT ALL DATA FILES
# ==============================================================================
cat("LOG [13/13]: Exporting Data Files...\n")

write.csv(t(final_scores), paste0(opt$out, "_Scores.csv"))
write.csv(t(z_res), paste0(opt$out, "_ZScores.csv"))
write.csv(t(gsva_res), paste0(opt$out, "_GSVA_Scores.csv"))
write.csv(meta, paste0(opt$out, "_Metadata.csv"))
write.csv(res_weighted$weights, paste0(opt$out, "_Weights_Weighted.csv"), row.names=FALSE)
write.csv(res_unweighted$weights, paste0(opt$out, "_Weights_Unweighted.csv"), row.names=FALSE)

writeLines(capture.output(sessionInfo()), paste0(dirname(opt$out), "/sessionInfo.txt"))

cat("\n╔═══════════════════════════════════════════════╗\n")
cat("║     ANALYSIS COMPLETE - v16.4 ULTIMATE        ║\n")
cat("╚═══════════════════════════════════════════════╝\n\n")
cat("✓ v16.2 Enhanced LLM reporting INTEGRATED\n")
cat("✓ v16.3 Comprehensive 2×4×N matrix INTEGRATED\n")
cat("✓ PCA gene drivers, plasticity stats, correlations\n")
cat("✓ Professional significance heatmap generated\n")
cat(sprintf("✓ Total tests per signature: %d\n", 2*4*(length(pairwise_combos)+1)))
cat(sprintf("✓ Trajectory patterns: %s\n", paste(names(pattern_table), "=", pattern_table, collapse=", ")))
cat(sprintf("✓ Consensus trends: %s\n", paste(names(consensus_table), "=", consensus_table, collapse=", ")))
cat(sprintf("\n📊 Significance matrix: %s\n", basename(paste0(opt$out, "_Significance_Matrix_Heatmap_mqc.png"))))
cat(sprintf("📝 Ultimate LLM summary: %s\n\n", basename(llm_file)))
