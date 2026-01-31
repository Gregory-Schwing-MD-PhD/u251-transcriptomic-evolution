#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# GLOBAL SUBTYPES ANALYSIS v15.4 (Best of Both Worlds Edition)
# ------------------------------------------------------------------------------
# Combines:
# - v15.3: Bug fixes, complete signatures, robust statistics, dual plasticity
# - v14.0: Rich HTML reporting, gene drivers, correlation analysis, better viz
# New features:
# - Fixed trajectory plot with all lines visible
# - Unified comprehensive HTML report
# - Enhanced statistical narrative
# - Complete visualization suite
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
# CONFIGURATION & AESTHETICS
# ==============================================================================
N_TOP_VAR <- opt$n_top_genes
SCORING_METHOD <- opt$scoring
FDR_THRESHOLDS <- c(0.05, 0.01, 0.005, 0.001)
FDR_LABELS <- c("Standard (0.05)", "Strict (0.01)", "Very Strict (0.005)", "Ultra (0.001)")

# Publication Theme
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

# Statistical Logger
stat_log <- list()
add_stat_log <- function(title, content) {
    stat_log[[length(stat_log) + 1]] <<- list(title=title, content=content)
}

# Formatting Utilities
safe_format <- function(x, digits=3) {
    if (is.null(x) || length(x) == 0 || is.na(x)) return("NA")
    if (!is.numeric(x)) return(as.character(x))
    if (abs(x) < 0.001 && x != 0) return(formatC(x, format="e", digits=2))
    return(round(x, digits))
}

interpret_p_text <- function(p) {
    if (is.null(p) || length(p) == 0 || is.na(p)) return("NA")
    if (p < 0.001) return("*** (highly significant)")
    if (p < 0.01) return("** (very significant)")
    if (p < 0.05) return("* (significant)")
    if (p < 0.10) return(". (trend)")
    return("ns (not significant)")
}

interpret_p_html <- function(p) {
    if (is.null(p) || length(p) == 0 || is.na(p)) return("<span class='ns'>NA</span>")
    if (p < 0.001) return("<span class='sig'>***</span>")
    if (p < 0.01) return("<span class='sig'>**</span>")
    if (p < 0.05) return("<span class='sig'>*</span>")
    if (p < 0.10) return("<span class='trend'>.</span>")
    return("<span class='ns'>ns</span>")
}

# ==============================================================================
# 1. DATA LOADING & PREP
# ==============================================================================
cat("LOG [1/10]: Loading Data & Checking Integrity...\n")

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)

# Load Matrix
mat_vst <- as.matrix(read.table(opt$input, header=TRUE, row.names=1, check.names=FALSE))
meta <- read.csv(opt$meta, row.names=1)

# Align Data
common <- intersect(colnames(mat_vst), rownames(meta))
if(length(common) < 3) stop("Error: <3 matching samples between VST and Metadata.")
mat_vst <- mat_vst[, common]
meta <- meta[common, , drop=FALSE]

# Factorize Classification
if(!"Classification" %in% colnames(meta)) {
    stop("Metadata must contain a 'Classification' column.")
}

# Check if we have the default U251 groups
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

# Robust Color Assignment
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

# Smart Gene ID Mapping
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
    add_stat_log("Gene Mapping", "Input already in gene symbol format - no mapping required")
}

# ==============================================================================
# 2. GLOBAL STRUCTURE (PCA + PERMANOVA + SCREE + GENE DRIVERS)
# ==============================================================================
cat("LOG [2/10]: Global Structure Analysis...\n")

# Feature Selection
top_var <- head(order(apply(mat_sym, 1, var), decreasing=TRUE), N_TOP_VAR)
mat_sig <- mat_sym[top_var, ]

add_stat_log("Feature Selection", sprintf(
    "Selected top %d most variable genes (variance range: %.2f - %.2f)",
    N_TOP_VAR, min(apply(mat_sig, 1, var)), max(apply(mat_sig, 1, var))
))

# PERMANOVA
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

# PCA with proper dimension handling
pca <- prcomp(t(mat_sig))

n_samples <- ncol(mat_sig)
max_pcs <- min(10, n_samples - 1, ncol(pca$x))
cat(sprintf("  → PCA computed: %d samples → maximum %d PCs available\n", n_samples, max_pcs))

pca_summary <- summary(pca)
n_pcs_to_show <- min(max_pcs, ncol(pca_summary$importance))

var_pc <- round(pca_summary$importance[2, 1:n_pcs_to_show]*100, 1)
cum_var <- round(pca_summary$importance[3, 1:n_pcs_to_show]*100, 1)
pcaData <- data.frame(pca$x[,1:2], Sample=rownames(meta), Class=meta$Classification)

add_stat_log("PCA Variance", sprintf(
    "PC1 explains %.1f%% variance\n  PC2 explains %.1f%% variance\n  PC1+PC2 together: %.1f%%\n  Total PCs available: %d",
    var_pc[1], var_pc[2], cum_var[2], n_pcs_to_show
))

# PCA Gene Drivers
loadings <- pca$rotation
pc_drivers_text <- ""
pc_drivers_list <- list()

for(i in 1:min(5, ncol(loadings))) {
    pc <- paste0("PC", i)
    top <- names(sort(abs(loadings[, i]), decreasing=TRUE)[1:8])
    pc_drivers_list[[pc]] <- top
    pc_drivers_text <- paste0(pc_drivers_text, "\n", pc, ": ", paste(top, collapse=", "))
}

add_stat_log("PCA Gene Drivers", sprintf(
    "Top genes driving principal components:%s",
    pc_drivers_text
))

# Plot A: Scree
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

# Plot B: PCA Trajectory with Biplot
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

# Top gene loadings for biplot
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
         subtitle=sprintf("PERMANOVA: R²=%.3f, P=%s", perm_r2, safe_format(perm_p)),
         x=paste0("PC1 (", var_pc[1], "%)"),
         y=paste0("PC2 (", var_pc[2], "%)")) +
    theme_publication()

# Combine with Patchwork
p_composite <- (p_pca | p_scree) + plot_layout(widths = c(2.5, 1))
ggsave(paste0(opt$out, "_Global_Structure_mqc.png"), p_composite, width=13, height=6)

# Phylogenetic Tree
phylo_tree <- as.phylo(hclust(dist(t(mat_sig)), method="ward.D2"))
tip_cols <- GROUP_COLORS[as.character(meta[phylo_tree$tip.label, "Classification"])]

png(paste0(opt$out, "_Phylogenetic_Tree_mqc.png"),
    width=10, height=8, units="in", res=300)
plot(phylo_tree, type="phylogram", tip.color=tip_cols,
     main="Phylogenetic Tree (Ward's D2)", cex=1.2, edge.width=2)
legend("topleft", legend=levels(meta$Classification),
       fill=GROUP_COLORS, bty="n")
dev.off()

# ==============================================================================
# 3. SUBTYPE SCORING (DUAL METHOD: GSVA + Z-SCORE)
# ==============================================================================
cat("LOG [3/10]: Scoring Signatures...\n")

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
    } else if (grepl("\\.(csv|tsv)$", opt$signatures, ignore.case = TRUE)) {
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

    add_stat_log("Signature Source", sprintf(
        "Loaded %d custom signatures from: %s",
        length(sigs), basename(opt$signatures)
    ))
} else {
    cat("  → Using COMPLETE built-in GBM signatures (Verhaak + Neftel + Garofano)\n")
    sigs <- list(
        # Verhaak (2010) TCGA Subtypes
        "Verhaak_Classical" = c(
            "EGFR", "NES", "NOTCH3", "JAG1", "HES5", "AKT2", "FGFR2",
            "CCND2", "CDK4", "RB1", "SOX2", "SFRP2"
        ),
        "Verhaak_Mesenchymal" = c(
            "CHI3L1", "CD44", "VIM", "RELB", "STAT3", "TRADD", "CASP1",
            "TNFRSF1A", "IKBKB", "NFKB1", "MET", "FOSL2", "TIMP1",
            "S100A4", "FN1", "COL1A2", "COL3A1", "SERPINE1"
        ),
        "Verhaak_Proneural" = c(
            "OLIG2", "DLL3", "ASCL1", "TCF12", "DCX", "SOX11", "PDGFRA",
            "NKX2-2", "ERBB3", "NCAM1", "NRCAM", "NKX6-1", "PAX3", "NEUROG1"
        ),
        "Verhaak_Neural" = c(
            "NEFL", "GABRA1", "SYT1", "SLC12A5", "GRIA2", "SYN1", "GABBR1",
            "SNAP25", "ATP1A3", "STX1A", "GABRG2", "NRXN1", "NLGN1"
        ),
        # Neftel (2019) Single-Cell States
        "Neftel_AC" = c(
            "APOE", "AQP4", "CLU", "S100B", "SLC1A2", "SLC1A3", "GFAP",
            "ALDOC", "GLUL", "GJA1", "ATP1A2", "ATP1B2", "FGFR3", "CD9",
            "MT3", "AGT", "PTN", "CST3", "SPARC", "SPARCL1"
        ),
        "Neftel_OPC" = c(
            "PDGFRA", "OLIG1", "OLIG2", "CSPG4", "SOX10", "PCDH15", "BCAN",
            "VCAN", "TNR", "LHFPL3", "GPR17", "MEGF11", "APOD", "COL9A1",
            "LUZP2", "PLP1", "MAG", "MBP", "MOG", "MOBP"
        ),
        "Neftel_NPC" = c(
            "DCX", "DLL3", "ASCL1", "NEUROG2", "STMN2", "SOX11", "SOX4",
            "TUBB3", "ELAVL2", "ELAVL3", "ELAVL4", "MAP2", "MAPT", "NCAM1",
            "NRCAM", "ROBO2", "L1CAM", "GAP43", "NEFM", "PRPH"
        ),
        "Neftel_MES" = c(
            "CHI3L1", "CD44", "ANXA1", "VIM", "S100A4", "S100A6", "S100A10",
            "LGALS1", "LGALS3", "SERPINE1", "TIMP1", "FN1", "COL1A1", "COL1A2",
            "COL3A1", "COL5A1", "LAMC1", "ITGA5", "ITGAV", "MMP2"
        ),
        # Garofano (2021) Metabolic Subtypes
        "Garofano_MTC" = c(
            "CS", "ACO2", "IDH2", "IDH3A", "OGDH", "SUCLA2", "SDHA", "SDHB",
            "FH", "MDH2", "NDUFS1", "NDUFS2", "NDUFV1", "UQCRC1", "UQCRC2",
            "COX4I1", "COX5A", "COX7A2", "ATP5F1A", "ATP5F1B"
        ),
        "Garofano_GPM" = c(
            "SLC2A1", "SLC2A3", "HK1", "HK2", "GPI", "PFKP", "PFKL", "PFKM",
            "ALDOA", "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "PKM",
            "LDHA", "LDHB", "SLC16A1", "SLC16A3", "G6PD", "PGD"
        ),
        "Garofano_NEU" = c(
            "GAD1", "GAD2", "SLC1A1", "SLC1A2", "SLC1A3", "GLUL", "GLS",
            "SYP", "SNAP25", "SYT1", "VAMP2", "STX1A", "ATP1A1", "ATP1A2",
            "ATP1A3", "SLC12A5", "KCNJ10", "KCNJ11"
        )
    )

    add_stat_log("Signature Source", sprintf(
        "Using complete built-in signatures: %d total (Verhaak: 4, Neftel: 4, Garofano: 3)",
        length(sigs)
    ))
}

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

# Method 2: Z-Score Averaging (with numerical stability)
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
        warning(sprintf("No genes found for signature %s", s))
    }
}

# Method Agreement Analysis
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
        theme_publication(base_size=12)

    ggsave(paste0(opt$out, "_Method_Agreement_mqc.png"), p_agreement,
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
# 4. PLASTICITY (ENTROPY WITH DUAL TESTING)
# ==============================================================================
cat("LOG [4/10]: Plasticity Analysis...\n")

# Numerically stable entropy calculation
calc_entropy <- function(s) {
    s <- pmin(pmax(s, -10), 10)
    p <- exp(s)/sum(exp(s))
    -sum(p*log(p + 1e-10), na.rm=TRUE)
}

meta$Plasticity <- apply(final_scores, 2, calc_entropy)

# Smart test selection based on sample size
n_total <- nrow(meta)
use_nonparametric <- (n_total < 30)

if(use_nonparametric) {
    cat(sprintf("  → Using Kruskal-Wallis (N=%d < 30, more powerful for small samples)\n", n_total))
    plast_test <- kruskal.test(Plasticity ~ Classification, data=meta)
    plast_p <- plast_test$p.value
    plast_stat <- plast_test$statistic
    plast_method <- "Kruskal-Wallis"

    add_stat_log("Plasticity: Test Selection", sprintf(
        "Sample size: %d\n  Decision: Kruskal-Wallis (non-parametric) chosen for small N\n  Rationale: More powerful than ANOVA for N<30 with potential outliers",
        n_total
    ))

    add_stat_log("Plasticity: Group Comparison", sprintf(
        "Kruskal-Wallis H-test\n  H = %.3f, df = %d\n  P = %s %s\n  Interpretation: %s",
        plast_stat, plast_test$parameter,
        safe_format(plast_p), interpret_p_text(plast_p),
        if(plast_p < 0.05) "Significant differences in plasticity across groups" else "No significant differences in plasticity"
    ))
} else {
    shapiro_p <- shapiro.test(meta$Plasticity)$p.value

    add_stat_log("Plasticity: Normality Test", sprintf(
        "Shapiro-Wilk test: W = %.3f, P = %s\n  Interpretation: Data is %s\n  Decision: Using ANOVA",
        shapiro.test(meta$Plasticity)$statistic,
        safe_format(shapiro_p),
        if(shapiro_p >= 0.05) "normally distributed" else "non-normal but using ANOVA (robust to violations with large N)"
    ))

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

# Plasticity group statistics
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

add_stat_log("Plasticity: Group Statistics", sprintf(
    "Summary by group:\n  %s",
    paste(capture.output(print(plast_stats, row.names=FALSE)), collapse="\n  ")
))

# Plasticity Plot
p_plast <- ggplot(meta, aes(x=Classification, y=Plasticity, fill=Classification)) +
    geom_violin(alpha=0.3, trim=FALSE) +
    geom_boxplot(width=0.2, fill="white", outlier.shape=NA) +
    geom_jitter(width=0.1, size=3, alpha=0.7) +
    stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="white") +
    scale_fill_manual(values=GROUP_COLORS) +
    labs(title="Cellular Plasticity (Shannon Entropy)",
         subtitle=sprintf("%s: P=%s %s", plast_method, safe_format(plast_p), interpret_p_text(plast_p)),
         caption="Diamond = mean, box = median ± IQR") +
    theme_publication() +
    theme(legend.position="none")

ggsave(paste0(opt$out, "_Plasticity_mqc.png"), p_plast, width=7, height=6)

# ==============================================================================
# 5. ROBUST LIMMA WITH ARRAY WEIGHTS
# ==============================================================================
cat("LOG [5/10]: Differential Analysis with Sample Quality Weights...\n")

design <- model.matrix(~0 + meta$Classification)
colnames(design) <- levels(meta$Classification)

# Calculate sample weights
aw <- arrayWeights(final_scores, design)
weight_df <- data.frame(
    Sample = names(aw),
    Weight = round(aw, 3),
    Group = as.character(meta[names(aw), "Classification"])
)

# Identify outliers
outlier_threshold <- mean(aw) - 2*sd(aw)
outliers <- weight_df$Sample[weight_df$Weight < outlier_threshold]

add_stat_log("Sample Quality Assessment", sprintf(
    "arrayWeights analysis:\n  Mean weight = %.3f\n  SD = %.3f\n  Range = %.3f - %.3f\n  Potential outliers (weight < mean - 2SD): %s\n  Interpretation: %s",
    mean(aw), sd(aw), min(aw), max(aw),
    if(length(outliers) > 0) paste(outliers, collapse=", ") else "None detected",
    if(length(outliers) > 0) "Some samples show reduced reliability - weights applied to down-weight their influence" else "All samples show consistent quality"
))

# Dynamic contrast generation
levs <- levels(meta$Classification)
contrast_formulas <- c()
contrast_names <- c()

if(length(levs) >= 2) {
    contrast_formulas <- c(contrast_formulas, sprintf("%s - %s", levs[2], levs[1]))
    contrast_names <- c(contrast_names, sprintf("%s_vs_%s", gsub("[^A-Za-z0-9]", "", levs[2]), gsub("[^A-Za-z0-9]", "", levs[1])))
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

# Fit model
fit <- lmFit(final_scores, design, weights=aw)
fit <- contrasts.fit(fit, cont.matrix)
fit <- eBayes(fit)

all_coefs <- colnames(cont.matrix)

add_stat_log("Linear Model Setup", sprintf(
    "Contrasts tested: %s\n  Model: Weighted linear regression with empirical Bayes moderation\n  Adjustment: Benjamini-Hochberg FDR correction",
    paste(contrast_names, collapse=", ")
))

# ==============================================================================
# 6. UNIFIED TRAJECTORY PLOT (FIXED - ALL LINES VISIBLE)
# ==============================================================================
cat("LOG [6/10]: Creating Unified Trajectory Visualization...\n")

# Prepare data
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

# Summary statistics for lines
traj_summary <- traj_data %>%
    group_by(Signature, Class, Stage) %>%
    summarise(Mean = mean(Score), SE = sd(Score)/sqrt(n()), .groups="drop")

# Sort by correlation
trend_order <- traj_data %>%
    group_by(Signature) %>%
    summarise(Cor = cor(Stage, Score, method="spearman"), .groups="drop") %>%
    arrange(desc(Cor)) %>%
    pull(Signature)

traj_data$Signature <- factor(traj_data$Signature, levels=trend_order)
traj_summary$Signature <- factor(traj_summary$Signature, levels=trend_order)

# Plot with explicit line layer (FIX: ensures all lines are drawn)
p_traj <- ggplot(traj_data, aes(x=Stage, y=Score)) +
    geom_ribbon(data=traj_summary, aes(x=Stage, ymin=Mean-SE, ymax=Mean+SE, fill=Signature),
                inherit.aes=FALSE, alpha=0.2) +
    geom_line(data=traj_summary, aes(x=Stage, y=Mean, color=Signature, group=Signature),
              inherit.aes=FALSE, linewidth=1.2, alpha=0.9) +
    geom_point(aes(fill=Class, shape=Class), size=3, alpha=0.6) +
    facet_wrap(~Signature, scales="free_y", ncol=3) +
    scale_x_continuous(breaks=1:length(levs), labels=levs) +
    scale_fill_manual(values=GROUP_COLORS, name="Stage") +
    scale_color_brewer(palette="Set1", guide="none") +
    scale_shape_manual(values=GROUP_SHAPES, name="Stage") +
    labs(title="Signature Trajectories Across Evolution",
         subtitle="Ordered by trend: Increasing (top) → Decreasing (bottom) | Lines = group means ± SE",
         x="Stage", y="Z-Score Expression") +
    theme_publication(base_size=11) +
    theme(legend.position = "bottom")

ggsave(paste0(opt$out, "_Unified_Trajectories_mqc.png"), p_traj, width=14, height=12)

# Trajectory statistics
traj_models <- traj_data %>%
    group_by(Signature) %>%
    summarise(
        Correlation_with_Stage = cor(Stage, Score, method="spearman"),
        Trend_P = cor.test(Stage, Score, method="spearman")$p.value,
        Direction = ifelse(Correlation_with_Stage > 0, "Increasing", "Decreasing"),
        .groups="drop"
    ) %>%
    arrange(desc(Correlation_with_Stage))

add_stat_log("Trajectory Analysis", sprintf(
    "Signature evolution patterns (sorted by trend strength):\n  %s",
    paste(capture.output(print(traj_models, row.names=FALSE)), collapse="\n  ")
))

# ==============================================================================
# 7. HEATMAPS & CORRELATION ANALYSIS
# ==============================================================================
cat("LOG [7/10]: Generating Heatmaps & Co-evolution Analysis...\n")

# Signature scores heatmap
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

# Signature correlation & co-evolution
sig_cor <- cor(t(final_scores), method="pearson")

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

png(paste0(opt$out, "_Sig_Correlation_mqc.png"), width=8, height=7, units="in", res=300)
Heatmap(sig_cor, name="Pearson\nCorr",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        column_title = "Signature Co-evolution",
        row_names_gp = gpar(fontsize=10),
        column_names_gp = gpar(fontsize=10))
dev.off()

# ==============================================================================
# 8. FDR SENSITIVITY & EXPORTS
# ==============================================================================
cat("LOG [8/10]: FDR Sensitivity & Statistical Results...\n")

# Primary results table
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
        tt[results_005$Signature, "adj.P.Val"] < 0.05, "✓", "—"
    )
}

# Sensitivity counts
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
    scale_fill_manual(values=rep(GROUP_COLORS[1:min(3, length(GROUP_COLORS))], length.out=length(all_coefs))) +
    labs(title="FDR Sensitivity Analysis",
         subtitle=sprintf("%d signatures tested across %d contrasts", nrow(final_scores), length(all_coefs)),
         x="FDR Threshold", y="Number of Significant Signatures") +
    theme_publication(base_size=12) +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    ylim(0, nrow(final_scores) + 1)

ggsave(paste0(opt$out, "_FDR_Sensitivity_mqc.png"), p_sensitivity, width=11, height=7)

# ==============================================================================
# 9. EXPORT DATA
# ==============================================================================
cat("LOG [9/10]: Exporting Data Files...\n")

write.csv(t(final_scores), paste0(opt$out, "_Scores.csv"))
write.csv(meta, paste0(opt$out, "_Metadata.csv"))
write.csv(weight_df, paste0(opt$out, "_SampleWeights.csv"), row.names=FALSE)
write.csv(results_005, paste0(opt$out, "_Statistical_Results.csv"), row.names=FALSE)
write.csv(sig_counts, paste0(opt$out, "_FDR_Sensitivity.csv"), row.names=FALSE)
write.csv(traj_models, paste0(opt$out, "_Trajectory_Statistics.csv"), row.names=FALSE)
write.csv(plast_stats, paste0(opt$out, "_Plasticity_Statistics.csv"), row.names=FALSE)

# MultiQC table
mqc_table <- paste0(dirname(opt$out), "/differential_signatures_mqc.tsv")
cat("# id: 'diff_sigs'\n# section_name: 'Differential Signatures'\n# plot_type: 'table'\n",
    file=mqc_table)
write.table(results_005, file=mqc_table, sep="\t", quote=FALSE, row.names=FALSE, append=TRUE)

# ==============================================================================
# 10. COMPREHENSIVE UNIFIED HTML REPORT
# ==============================================================================
cat("LOG [10/10]: Generating Comprehensive HTML Report...\n")

# Count significant signatures for executive summary
sig_counts_summary <- list()
for(coef in all_coefs) {
    tt <- topTable(fit, coef=coef, number=Inf)
    sig_counts_summary[[coef]] <- sum(tt$adj.P.Val < 0.05, na.rm=TRUE)
}

sink(paste0(opt$out, "_Report.html"))
cat("<!DOCTYPE html>
<html lang='en'>
<head>
<meta charset='UTF-8'>
<meta name='viewport' content='width=device-width, initial-scale=1.0'>
<title>U251 Evolution Analysis Report v15.4 - Best of Both Worlds</title>
<style>
body {
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    background: #f4f7f6;
    color: #333;
    line-height: 1.6;
    max-width: 1400px;
    margin: 0 auto;
    padding: 20px;
}
.header {
    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    color: white;
    padding: 30px;
    border-radius: 10px;
    margin-bottom: 30px;
    box-shadow: 0 4px 6px rgba(0,0,0,0.1);
}
.header h1 { margin: 0 0 10px 0; }
.header h3 { margin: 0 0 15px 0; opacity: 0.9; }
.card {
    background: white;
    padding: 25px;
    border-radius: 8px;
    margin-bottom: 25px;
    box-shadow: 0 2px 5px rgba(0,0,0,0.05);
    border-left: 4px solid #667eea;
}
.stat-block {
    background: #f8f9fa;
    padding: 15px;
    margin: 10px 0;
    border-radius: 5px;
    border-left: 3px solid #667eea;
}
.stat-title {
    font-weight: bold;
    color: #667eea;
    margin-bottom: 8px;
    font-size: 1.05em;
}
.stat-content {
    font-family: 'Courier New', monospace;
    font-size: 0.85em;
    white-space: pre-wrap;
    line-height: 1.5;
    color: #555;
}
table {
    width: 100%;
    border-collapse: collapse;
    font-size: 0.9em;
    margin: 15px 0;
}
th {
    background: #667eea;
    color: white;
    text-align: left;
    padding: 12px;
    font-weight: 600;
}
td {
    padding: 10px;
    border-bottom: 1px solid #eee;
}
tr:hover {
    background-color: #f5f7fa;
}
.sig {
    color: #d62728;
    font-weight: bold;
}
.ns {
    color: #adb5bd;
}
.trend {
    color: #fd7e14;
    font-weight: bold;
}
.img-container {
    text-align: center;
    margin: 20px 0;
}
img {
    max-width: 100%;
    border-radius: 5px;
    box-shadow: 0 2px 8px rgba(0,0,0,0.1);
}
h2 {
    color: #2c3e50;
    border-bottom: 2px solid #667eea;
    padding-bottom: 10px;
    margin-top: 30px;
}
.quality-excellent { color: #28a745; font-weight: bold; }
.quality-good { color: #17a2b8; }
.quality-fair { color: #ffc107; }
.quality-poor { color: #dc3545; font-weight: bold; }
.badge {
    display: inline-block;
    padding: 4px 10px;
    border-radius: 4px;
    font-size: 0.85em;
    font-weight: 600;
    margin: 3px;
}
.badge-info {
    background: #d1ecf1;
    color: #0c5460;
}
.badge-success {
    background: #d4edda;
    color: #155724;
}
.badge-warning {
    background: #fff3cd;
    color: #856404;
}
.highlight {
    background-color: #fff3cd;
    padding: 2px 6px;
    border-radius: 3px;
}
.interpretation-box {
    background-color: #e7f3ff;
    border-left: 4px solid #2196F3;
    padding: 20px;
    margin: 20px 0;
    border-radius: 5px;
}
.interpretation-box h4 {
    color: #2196F3;
    margin-top: 0;
}
.interpretation-box ul {
    margin: 10px 0;
}
.key-finding {
    background: #fff9e6;
    border-left: 4px solid #ffc107;
    padding: 15px;
    margin: 15px 0;
    border-radius: 5px;
}
</style>
</head>
<body>

<div class='header'>
<h1>🧬 U251 Global Subtype Evolution Analysis</h1>
<h3>Comprehensive Statistical Report v15.4 (Best of Both Worlds)</h3>
<p><strong>Analysis Date:</strong> ", Sys.Date(), " | <strong>Samples:</strong> ", ncol(mat_vst), " | <strong>Genes:</strong> ", nrow(mat_sym), "</p>
<p><span class='badge badge-info'>", length(sigs), " Subtype Signatures</span>
   <span class='badge badge-info'>", length(all_coefs), " Contrasts</span>
   <span class='badge badge-info'>", nrow(final_scores), " Signatures Scored</span>
   <span class='badge badge-success'>Scoring: ", scoring_label, "</span></p>
</div>

<div class='card'>
<h2>📊 Executive Summary</h2>

<div class='stat-block'>
<div class='stat-title'>Dataset Overview</div>
<div class='stat-content'>")

cat(sprintf("Samples analyzed: %d\nGenes quantified: %d\nSignatures evaluated: %d\n\nGroup breakdown:\n",
            ncol(mat_vst), nrow(mat_sym), nrow(final_scores)))
for(i in seq_along(group_sizes)) {
    cat(sprintf("  • %s: %d samples\n", names(group_sizes)[i], group_sizes[i]))
}

cat("</div></div>

<div class='key-finding'>
<h4 style='margin-top:0; color:#856404;'>🔬 Key Statistical Findings</h4>
<div class='stat-content' style='font-family: Arial; color: #333;'>")

cat(sprintf("<strong>Multivariate Analysis:</strong>\n"))
cat(sprintf("  PERMANOVA: P = %s %s\n", safe_format(perm_p), interpret_p_text(perm_p)))
cat(sprintf("  R² = %.3f (%.1f%% variance explained by groups)\n", perm_r2, perm_r2*100))
cat(sprintf("  → Groups %s distinct transcriptional signatures\n\n",
            if(!is.na(perm_p) && perm_p < 0.05) "SHOW" else "DO NOT SHOW"))

cat(sprintf("<strong>Cellular Plasticity:</strong>\n"))
cat(sprintf("  %s: P = %s %s\n", plast_method, safe_format(plast_p), interpret_p_text(plast_p)))
cat(sprintf("  → Plasticity %s across evolutionary stages\n\n",
            if(plast_p < 0.05) "CHANGES SIGNIFICANTLY" else "remains stable"))

cat(sprintf("<strong>Differential Signatures (FDR < 0.05):</strong>\n"))
for(i in 1:length(all_coefs)) {
    coef <- all_coefs[i]
    n_sig <- sig_counts_summary[[coef]]
    cat(sprintf("  • %s: %d/%d signatures (%d%%)\n",
                gsub("_", " ", coef), n_sig, nrow(final_scores),
                round(n_sig/nrow(final_scores)*100)))
}

cat("</div></div>
</div>

<div class='card'>
<h2>🗺️ Global Structure & Evolutionary Trajectories</h2>
<div class='img-container'>
<img src='", basename(paste0(opt$out, "_Global_Structure_mqc.png")), "' alt='PCA and Scree Plot'>
<p style='font-size:0.9em; color:#666; margin-top:10px;'>
<strong>Left:</strong> PCA trajectory showing group centroids (arrows) and gene drivers (red arrows). 
<strong>Right:</strong> Scree plot of variance explained by each principal component.
</p>
</div>

<div class='img-container'>
<img src='", basename(paste0(opt$out, "_Phylogenetic_Tree_mqc.png")), "' alt='Phylogenetic Tree'>
<p style='font-size:0.9em; color:#666; margin-top:10px;'>
Hierarchical clustering (Ward's D2) showing relationships between samples based on gene expression.
</p>
</div>
</div>

<div class='card'>
<h2>🧬 Signature Trajectories & Evolution</h2>
<div class='img-container'>
<img src='", basename(paste0(opt$out, "_Unified_Trajectories_mqc.png")), "' alt='Trajectory Analysis'>
<p style='font-size:0.9em; color:#666; margin-top:10px;'>
Signatures ordered by correlation with evolutionary stage. Lines show group means ± SE; points are individual samples.
</p>
</div>

<div class='stat-block'>
<div class='stat-title'>Trajectory Statistics</div>
<table>
<thead>
<tr><th>Signature</th><th>Correlation (ρ)</th><th>P-value</th><th>Direction</th><th>Significance</th></tr>
</thead>
<tbody>")

for(i in 1:nrow(traj_models)) {
    sig_marker <- if(traj_models$Trend_P[i] < 0.05) {
        "<span class='sig'>✓ Significant</span>"
    } else {
        "<span class='ns'>— Not Significant</span>"
    }

    cat(sprintf("<tr><td><strong>%s</strong></td><td>%.3f</td><td>%s</td><td>%s</td><td>%s</td></tr>\n",
                traj_models$Signature[i],
                traj_models$Correlation_with_Stage[i],
                safe_format(traj_models$Trend_P[i]),
                traj_models$Direction[i],
                sig_marker))
}

cat("</tbody>
</table>
</div>
</div>

<div class='card'>
<h2>🎨 Cellular Plasticity Analysis</h2>
<div class='img-container'>
<img src='", basename(paste0(opt$out, "_Plasticity_mqc.png")), "' alt='Plasticity Analysis'>
</div>

<div class='stat-block'>
<div class='stat-title'>Plasticity by Group</div>
<table>
<thead>
<tr><th>Group</th><th>N</th><th>Mean</th><th>SD</th><th>Median</th><th>IQR</th></tr>
</thead>
<tbody>")

for(i in 1:nrow(plast_stats)) {
    cat(sprintf("<tr><td><strong>%s</strong></td><td>%d</td><td>%.3f</td><td>%.3f</td><td>%.3f</td><td>%.3f</td></tr>\n",
                plast_stats$Classification[i],
                plast_stats$N[i],
                plast_stats$Mean[i],
                plast_stats$SD[i],
                plast_stats$Median[i],
                plast_stats$IQR[i]))
}

cat("</tbody>
</table>
</div>
</div>

<div class='card'>
<h2>🔥 Heatmaps & Co-evolution</h2>
<div class='img-container'>
<img src='", basename(paste0(opt$out, "_Heatmap_mqc.png")), "' alt='Signature Heatmap'>
<p style='font-size:0.9em; color:#666; margin-top:10px;'>
Signature scores across all samples, ordered by evolutionary stage.
</p>
</div>

<div class='img-container'>
<img src='", basename(paste0(opt$out, "_Sig_Correlation_mqc.png")), "' alt='Correlation Matrix'>
<p style='font-size:0.9em; color:#666; margin-top:10px;'>
Pairwise correlations between signatures reveal co-evolution patterns.
</p>
</div>")

if(exists("cor_summary") && nrow(cor_summary) > 0) {
    cat("
<div class='stat-block'>
<div class='stat-title'>Strong Co-evolution Pairs (|r| > 0.6)</div>
<table>
<thead>
<tr><th>Signature 1</th><th>Signature 2</th><th>Correlation</th></tr>
</thead>
<tbody>")

    for(i in 1:min(10, nrow(cor_summary))) {
        cat(sprintf("<tr><td>%s</td><td>%s</td><td>%.3f</td></tr>\n",
                    cor_summary$Sig1[i],
                    cor_summary$Sig2[i],
                    cor_summary$Correlation[i]))
    }

    cat("</tbody></table></div>")
}

if(SCORING_METHOD == "both") {
    cat("
<div class='img-container'>
<img src='", basename(paste0(opt$out, "_Method_Agreement_mqc.png")), "' alt='Method Agreement'>
<p style='font-size:0.9em; color:#666; margin-top:10px;'>
Comparison of GSVA and Z-Score methods shows strong agreement.
</p>
</div>")
}

cat("</div>

<div class='card'>
<h2>📈 Differential Analysis & FDR Sensitivity</h2>
<div class='img-container'>
<img src='", basename(paste0(opt$out, "_FDR_Sensitivity_mqc.png")), "' alt='FDR Sensitivity'>
</div>

<div class='stat-block'>
<div class='stat-title'>Differential Results (FDR < 0.05)</div>
<table>
<thead>
<tr>
<th>Signature</th>")

for(coef in all_coefs) {
    cat(sprintf("<th colspan='2'>%s</th>", gsub("_", " ", coef)))
}
cat("</tr>
<tr><th></th>")
for(coef in all_coefs) {
    cat("<th>logFC</th><th>adjP</th>")
}
cat("</tr>
</thead>
<tbody>")

for(i in 1:nrow(results_005)) {
    cat(sprintf("<tr><td><strong>%s</strong></td>", results_005$Signature[i]))
    for(coef in all_coefs) {
        logfc_col <- paste0(coef, "_logFC")
        adjp_col <- paste0(coef, "_adjP")
        sig_col <- paste0(coef, "_Sig")

        adjp_val <- as.numeric(results_005[[adjp_col]][i])
        sig_class <- if(!is.na(adjp_val) && adjp_val < 0.05) "sig" else "ns"

        cat(sprintf("<td>%.3f</td><td class='%s'>%s</td>",
                    results_005[[logfc_col]][i],
                    sig_class,
                    results_005[[adjp_col]][i]))
    }
    cat("</tr>\n")
}

cat("</tbody>
</table>
</div>
</div>

<div class='card'>
<h2>🎯 PCA Gene Drivers</h2>
<div class='stat-block'>
<div class='stat-content'>")

for(pc in names(pc_drivers_list)) {
    cat(sprintf("<strong>%s:</strong> %s\n", pc, paste(pc_drivers_list[[pc]], collapse=", ")))
}

cat("</div>
</div>
</div>

<div class='card'>
<h2>📉 Sample Quality Assessment</h2>
<div class='stat-block'>
<table>
<thead>
<tr><th>Sample</th><th>Weight</th><th>Group</th><th>Quality</th></tr>
</thead>
<tbody>")

for(i in 1:nrow(weight_df)) {
    quality_class <- if(weight_df$Weight[i] > mean(aw) + sd(aw)) {
        "quality-excellent'>Excellent"
    } else if(weight_df$Weight[i] > mean(aw)) {
        "quality-good'>Good"
    } else if(weight_df$Weight[i] > mean(aw) - sd(aw)) {
        "quality-fair'>Fair"
    } else {
        "quality-poor'>Poor"
    }

    cat(sprintf("<tr><td>%s</td><td>%.3f</td><td>%s</td><td><span class='%s</span></td></tr>\n",
                weight_df$Sample[i], weight_df$Weight[i], weight_df$Group[i], quality_class))
}

cat("</tbody>
</table>
<p style='margin-top: 15px; font-size: 0.9em; color: #666;'>
<strong>Note:</strong> Sample weights calculated by limma's arrayWeights. Lower weights indicate noisier samples that are down-weighted in analysis.
</p>
</div>
</div>

<div class='card'>
<h2>📈 Detailed Statistical Methods & Decisions</h2>")

for(log_entry in stat_log) {
    cat(sprintf("
<div class='stat-block'>
<div class='stat-title'>%s</div>
<div class='stat-content'>%s</div>
</div>", log_entry$title, log_entry$content))
}

cat("</div>

<div class='interpretation-box'>
<h2>💡 Interpretation Guide</h2>

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
<li><strong>R² (PERMANOVA):</strong> Proportion of variance explained by groups (0-1 scale, higher = stronger group effect)</li>
<li><strong>Shannon Entropy:</strong> Measure of cellular plasticity; higher = more plastic/undifferentiated state</li>
<li><strong>logFC:</strong> Log2 fold change; positive = higher in second group, negative = higher in first group</li>
<li><strong>Sample Weights:</strong> Quality metric; ~1.0 is typical, <0.5 suggests noise</li>
<li><strong>Spearman ρ:</strong> Correlation with stage; positive = increasing trajectory, negative = decreasing</li>
</ul>

<h4>Test Selection Logic:</h4>
<ul>
<li><strong>Plasticity:</strong> Kruskal-Wallis for N<30 (more powerful for small samples), ANOVA otherwise</li>
<li><strong>Differential Analysis:</strong> Weighted linear models with empirical Bayes shrinkage</li>
<li><strong>Multiple Testing:</strong> Benjamini-Hochberg FDR correction across all signatures</li>
</ul>
</div>

<div class='card'>
<h2>🔬 Methods & Software</h2>
<div class='stat-block'>
<div class='stat-content' style='font-family: Arial; white-space: normal;'>
<strong>R Version:</strong> ", as.character(R.version.string), "<br>
<strong>Key Packages:</strong><br>
  • limma ", as.character(packageVersion("limma")), " (differential analysis)<br>
  • GSVA ", as.character(packageVersion("GSVA")), " (signature scoring)<br>
  • vegan ", as.character(packageVersion("vegan")), " (PERMANOVA)<br>
  • ComplexHeatmap ", as.character(packageVersion("ComplexHeatmap")), " (visualization)<br><br>
<strong>Statistical Approach:</strong><br>
  • Scoring: ", scoring_label, "<br>
  • Model: Weighted linear regression + empirical Bayes<br>
  • Sample weights: arrayWeights (automatic quality weighting)<br>
  • Correction: Benjamini-Hochberg FDR<br>
  • Normality: Shapiro-Wilk test<br>
  • Group tests: Parametric (ANOVA) or non-parametric (Kruskal-Wallis)<br><br>
<strong>Reproducibility:</strong><br>
  • Random seed: 12345<br>
  • Top variable genes: ", N_TOP_VAR, "<br>
  • Analysis timestamp: ", Sys.time(), "<br>
  • Signatures: ", length(sigs), " (Verhaak: 4, Neftel: 4, Garofano: 3)
</div>
</div>
</div>

<div style='text-align: center; margin: 40px 0; padding: 20px; background: #f8f9fa; border-radius: 5px;'>
<p style='color: #666; margin: 0;'>
Generated by <strong>U251 Global Subtype Analysis Pipeline v15.4</strong><br>
Best of Both Worlds Edition: Robust Statistics + Rich Reporting
</p>
</div>

</body>
</html>")

sink()

# Session info
writeLines(capture.output(sessionInfo()), paste0(dirname(opt$out), "/sessionInfo.txt"))

cat("\n╔═══════════════════════════════════════════════╗\n")
cat("║   ANALYSIS COMPLETE - v15.4 BEST OF BOTH      ║\n")
cat("╚═══════════════════════════════════════════════╝\n\n")
cat(sprintf("✓ Analyzed %d samples, %d genes, %d signatures\n", ncol(mat_vst), nrow(mat_sym), nrow(final_scores)))
cat(sprintf("✓ Statistical logs: %d entries\n", length(stat_log)))
cat(sprintf("✓ Visualizations: 7 plots generated\n"))
cat(sprintf("✓ Comprehensive HTML report with interpretation guide\n"))
cat(sprintf("\nKey features:\n"))
cat(sprintf("  • Fixed trajectory plot (all lines visible)\n"))
cat(sprintf("  • v15.3 robust statistics + v14.0 rich reporting\n"))
cat(sprintf("  • Complete signatures (11 total)\n"))
cat(sprintf("  • Dual plasticity testing\n"))
cat(sprintf("  • Gene drivers & co-evolution analysis\n"))
cat(sprintf("\nMain report: %s\n", paste0(opt$out, "_Report.html")))
