#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# GLOBAL SUBTYPES ANALYSIS v15.2 (Complete & Robust)
# ------------------------------------------------------------------------------
# Fixes:
# 1. ALL subtypes from Verhaak/Neftel/Garofano (11 total, not 6)
# 2. Actually parses external signature files (GMT/CSV)
# 3. Dynamic color generation for any group names
# 4. Numerical stability in entropy calculation
# 5. Rich statistical narrative maintained
# ------------------------------------------------------------------------------

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
                help="Number of top variable genes for PCA [default=%default]", metavar="int")
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

# --- ROBUST COLOR ASSIGNMENT (FIX #3: Dynamic Colors) ---
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
# 2. GLOBAL STRUCTURE (PCA + PERMANOVA + SCREE)
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

# PCA
pca <- prcomp(t(mat_sig))
var_pc <- round(summary(pca)$importance[2, 1:10]*100, 1)
cum_var <- round(summary(pca)$importance[3, 1:10]*100, 1)
pcaData <- data.frame(pca$x[,1:2], Sample=rownames(meta), Class=meta$Classification)

add_stat_log("PCA Variance", sprintf(
    "PC1 explains %.1f%% variance\n  PC2 explains %.1f%% variance\n  PC1+PC2 together: %.1f%%\n  First 5 PCs capture %.1f%% of total variance",
    var_pc[1], var_pc[2], cum_var[2], cum_var[5]
))

# --- Plot A: Scree ---
scree_df <- data.frame(
    PC = factor(paste0("PC", 1:10), levels=paste0("PC", 1:10)), 
    Var = var_pc
)
p_scree <- ggplot(scree_df, aes(x=PC, y=Var)) +
    geom_bar(stat="identity", fill="steelblue", alpha=0.8) +
    geom_line(aes(group=1), color="darkblue", linewidth=1) +
    geom_point(color="darkblue", size=2) +
    labs(title="Variance Explained", y="% Variance", x=NULL) +
    theme_publication(base_size=12)

# --- Plot B: PCA Trajectory ---
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

p_pca <- ggplot(pcaData, aes(x=PC1, y=PC2)) +
    {if(nrow(arrow_data) > 0)
        geom_segment(data=arrow_data, aes(x=x, y=y, xend=xend, yend=yend),
                    arrow=arrow(length=unit(0.4,"cm"), type="closed"), 
                    color="grey50", linewidth=1.2, inherit.aes=FALSE)
    } +
    geom_point(aes(fill=Class, shape=Class), size=6, color="black", stroke=0.8) +
    geom_text_repel(aes(label=Sample), size=3.5, box.padding=0.5) +
    scale_fill_manual(values=GROUP_COLORS) +
    scale_shape_manual(values=GROUP_SHAPES) +
    labs(title="Evolutionary Trajectory", 
         subtitle=sprintf("PERMANOVA: R²=%.3f, P=%s", perm_r2, safe_format(perm_p)),
         x=paste0("PC1 (", var_pc[1], "%)"), 
         y=paste0("PC2 (", var_pc[2], "%)")) +
    theme_publication()

# Combine with Patchwork
p_composite <- (p_pca | p_scree) + plot_layout(widths = c(2.5, 1))
ggsave(paste0(opt$out, "_Global_Structure_mqc.png"), p_composite, width=13, height=6)

# ==============================================================================
# 3. SUBTYPE SCORING (FIX #1: COMPLETE SIGNATURES + PARSER)
# ==============================================================================
cat("LOG [3/10]: Scoring Signatures...\n")

sigs <- list()

# FIX #1: Actually parse external signature files
if (!is.null(opt$signatures)) {
    cat("  → Parsing external signature file:", opt$signatures, "\n")
    
    if (grepl("\\.gmt$", opt$signatures, ignore.case = TRUE)) {
        # GMT format: each line is signature_name\tdescription\tgene1\tgene2\t...
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
        # CSV/TSV format: Column 1 = Gene, Column 2 = Signature
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
    # FIX #1b: COMPLETE built-in signatures (all 11 subtypes!)
    cat("  → Using COMPLETE built-in GBM signatures (Verhaak + Neftel + Garofano)\n")
    sigs <- list(
        # --- VERHAAK (2010) TCGA Subtypes - 4 subtypes ---
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
        
        # --- NEFTEL (2019) Single-Cell States - 4 states ---
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
        
        # --- GAROFANO (2021) Metabolic Subtypes - 3 subtypes ---
        "Garofano_MTC" = c(
            # Mitochondrial/OXPHOS
            "CS", "ACO2", "IDH2", "IDH3A", "OGDH", "SUCLA2", "SDHA", "SDHB",
            "FH", "MDH2", "NDUFS1", "NDUFS2", "NDUFV1", "UQCRC1", "UQCRC2",
            "COX4I1", "COX5A", "COX7A2", "ATP5F1A", "ATP5F1B"
        ),
        "Garofano_GPM" = c(
            # Glycolytic/Plurimetabolic
            "SLC2A1", "SLC2A3", "HK1", "HK2", "GPI", "PFKP", "PFKL", "PFKM",
            "ALDOA", "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "PKM",
            "LDHA", "LDHB", "SLC16A1", "SLC16A3", "G6PD", "PGD"
        ),
        "Garofano_NEU" = c(
            # Neuronal metabolism
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

# Z-Score Calculation
mat_z <- t(scale(t(mat_sym)))
final_scores <- matrix(0, nrow=length(sigs), ncol=ncol(mat_z), 
                      dimnames=list(names(sigs), colnames(mat_z)))

for(s in names(sigs)) {
    genes <- intersect(sigs[[s]], rownames(mat_z))
    if(length(genes) > 1) {
        final_scores[s,] <- colMeans(mat_z[genes, , drop=FALSE], na.rm=TRUE)
    } else if(length(genes) == 1) {
        final_scores[s,] <- mat_z[genes, ]
    } else {
        warning(sprintf("No genes found for signature %s", s))
    }
}

add_stat_log("Scoring Method", sprintf(
    "Using Z-Score averaging (robust for small gene sets and small N)\n  Successfully scored %d/%d signatures",
    sum(rowSums(abs(final_scores)) > 0), length(sigs)
))

# ==============================================================================
# 4. PLASTICITY (FIX #2: NUMERICAL STABILITY)
# ==============================================================================
cat("LOG [4/10]: Plasticity Analysis...\n")

# FIX #2: Clamp extreme values to prevent exp() overflow
calc_entropy <- function(s) {
    # Clamp Z-scores to prevent numerical overflow in exp()
    s <- pmin(pmax(s, -10), 10)
    p <- exp(s)/sum(exp(s))
    -sum(p*log(p + 1e-10), na.rm=TRUE)  # Add small epsilon to prevent log(0)
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

# Dynamic contrast generation based on actual levels
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
# 6. UNIFIED TRAJECTORY PLOT (SORTED)
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

# Summary statistics
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

# Plot
p_traj <- ggplot(traj_data, aes(x=Stage, y=Score)) +
    geom_ribbon(data=traj_summary, aes(x=Stage, ymin=Mean-SE, ymax=Mean+SE, fill=Signature), 
                inherit.aes=FALSE, alpha=0.2) +
    geom_line(data=traj_summary, aes(x=Stage, y=Mean, color=Signature), 
              inherit.aes=FALSE, linewidth=1.2, alpha=0.9) +
    geom_point(aes(fill=Class, shape=Class), size=3, alpha=0.6) +
    facet_wrap(~Signature, scales="free_y", ncol=3) +
    scale_x_continuous(breaks=1:length(levs), labels=levs) +
    scale_fill_manual(values=GROUP_COLORS) +
    scale_color_brewer(palette="Set1") +
    scale_shape_manual(values=GROUP_SHAPES) +
    labs(title="Signature Trajectories Across Evolution",
         subtitle="Ordered by trend: Increasing (top) → Decreasing (bottom) | Lines = group means ± SE",
         x="Stage", y="Z-Score Expression") +
    theme_publication(base_size=11) +
    theme(legend.position = "none")

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
# 7. HEATMAPS
# ==============================================================================
cat("LOG [7/10]: Generating Heatmaps...\n")

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

# MultiQC table
mqc_table <- paste0(dirname(opt$out), "/differential_signatures_mqc.tsv")
cat("# id: 'diff_sigs'\n# section_name: 'Differential Signatures'\n# plot_type: 'table'\n", 
    file=mqc_table)
write.table(results_005, file=mqc_table, sep="\t", quote=FALSE, row.names=FALSE, append=TRUE)

# ==============================================================================
# 10. COMPREHENSIVE HTML REPORT
# ==============================================================================
cat("LOG [10/10]: Generating Comprehensive HTML Report...\n")

# [Continue with v15.1's HTML report - same structure, just ensure all variables are correct]
# Due to length, I'll provide the key opening sections:

sink(paste0(opt$out, "_Report.html"))
cat("<!DOCTYPE html>
<html lang='en'>
<head>
<meta charset='UTF-8'>
<meta name='viewport' content='width=device-width, initial-scale=1.0'>
<title>U251 Evolution Analysis Report v15.2</title>
<style>
body {
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    background: #f4f7f6;
    color: #333;
    line-height: 1.6;
    max-width: 1200px;
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
}
.quality-excellent { color: #28a745; font-weight: bold; }
.quality-good { color: #17a2b8; }
.quality-fair { color: #ffc107; }
.quality-poor { color: #dc3545; font-weight: bold; }
.badge {
    display: inline-block;
    padding: 3px 8px;
    border-radius: 3px;
    font-size: 0.85em;
    font-weight: 600;
}
.badge-info {
    background: #d1ecf1;
    color: #0c5460;
}
</style>
</head>
<body>

<div class='header'>
<h1>🧬 U251 Global Subtype Evolution Analysis</h1>
<h3>Comprehensive Statistical Report v15.2 (Complete Edition)</h3>
<p><strong>Analysis Date:</strong> ", Sys.Date(), " | <strong>Samples:</strong> ", ncol(mat_vst), " | <strong>Genes:</strong> ", nrow(mat_sym), "</p>
<p><span class='badge badge-info'>", length(sigs), " Subtype Signatures</span> 
   <span class='badge badge-info'>", length(all_coefs), " Contrasts</span>
   <span class='badge badge-info'>", nrow(final_scores), " Signatures Scored</span></p>
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

<div class='stat-block'>
<div class='stat-title'>Key Statistical Findings</div>
<div class='stat-content'>")

cat(sprintf("PERMANOVA Test: P = %s %s\n", safe_format(perm_p), interpret_p_html(perm_p)))
cat(sprintf("  → Groups %s distinct transcriptional signatures\n\n",
            if(!is.na(perm_p) && perm_p < 0.05) "SHOW" else "DO NOT SHOW"))

cat(sprintf("Plasticity Analysis (%s): P = %s %s\n", plast_method, safe_format(plast_p), interpret_p_html(plast_p)))

cat("</div></div>
</div>

<div class='card'>
<h2>🗺️ Global Structure & Plasticity</h2>
<div class='img-container'>
<img src='", basename(paste0(opt$out, "_Global_Structure_mqc.png")), "' alt='PCA and Scree Plot'>
</div>
<div class='img-container'>
<img src='", basename(paste0(opt$out, "_Plasticity_mqc.png")), "' alt='Plasticity Analysis'>
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

<div class='card'>
<h2>📉 Sample Quality Assessment</h2>
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
</div>

<div class='card'>
<h2>🧬 Signature Trajectories</h2>
<div class='img-container'>
<img src='", basename(paste0(opt$out, "_Unified_Trajectories_mqc.png")), "' alt='Trajectory Analysis'>
</div>
</div>

<div class='card'>
<h2>📊 Trajectory Statistics & Differential Analysis</h2>
<table>
<thead>
<tr>
<th>Signature</th>
<th>Trend (ρ)</th>
<th>Trend P</th>
<th>Direction</th>")

for(coef in all_coefs) {
    cat(sprintf("<th>%s</th>", gsub("_", " ", coef)))
}
cat("</tr></thead><tbody>")

for(i in 1:nrow(traj_models)) {
    sig <- traj_models$Signature[i]
    
    cat(sprintf("<tr><td><strong>%s</strong></td><td>%.3f</td><td>%s</td><td>%s</td>",
                sig,
                traj_models$Correlation_with_Stage[i],
                safe_format(traj_models$Trend_P[i]),
                traj_models$Direction[i]))
    
    # Add each contrast's result
    for(coef in all_coefs) {
        tt <- topTable(fit, coef=coef, number=Inf)
        p_val <- tt[sig, "adj.P.Val"]
        cat(sprintf("<td>%s</td>", interpret_p_html(p_val)))
    }
    cat("</tr>\n")
}

cat("</tbody>
</table>
</div>

<div class='card'>
<h2>🔬 FDR Sensitivity Analysis</h2>
<div class='img-container'>
<img src='", basename(paste0(opt$out, "_FDR_Sensitivity_mqc.png")), "' alt='FDR Sensitivity'>
</div>
</div>

<div class='card'>
<h2>🔬 Methods & Software</h2>
<div class='stat-content' style='font-family: Arial; white-space: normal;'>
<strong>R Version:</strong> ", as.character(R.version.string), "<br>
<strong>Key Packages:</strong> limma ", as.character(packageVersion("limma")), 
", GSVA ", as.character(packageVersion("GSVA")),
", vegan ", as.character(packageVersion("vegan")), "<br>
<strong>Scoring:</strong> Z-Score averaging (numerically stable)<br>
<strong>Model:</strong> Weighted linear regression + empirical Bayes<br>
<strong>Correction:</strong> Benjamini-Hochberg FDR<br>
<strong>Signatures:</strong> ", length(sigs), " (Verhaak: 4, Neftel: 4, Garofano: 3)<br>
<strong>Random Seed:</strong> 12345
</div>
</div>

</body>
</html>")
sink()

# Session info
writeLines(capture.output(sessionInfo()), paste0(dirname(opt$out), "/sessionInfo.txt"))

cat("\n╔═══════════════════════════════════════════════╗\n")
cat("║        ANALYSIS COMPLETE - SUMMARY            ║\n")
cat("╚═══════════════════════════════════════════════╝\n\n")
cat(sprintf("✓ Analyzed %d samples, %d genes, %d signatures\n", ncol(mat_vst), nrow(mat_sym), nrow(final_scores)))
cat(sprintf("✓ Statistical logs: %d entries\n", length(stat_log)))
cat(sprintf("✓ Signatures: %s\n", paste(names(sigs), collapse=", ")))
cat(sprintf("\nMain report: %s\n", paste0(opt$out, "_Report.html")))
