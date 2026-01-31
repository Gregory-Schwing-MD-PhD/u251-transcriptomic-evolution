#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# GLOBAL SUBTYPES ANALYSIS v15.5 (ULTIMATE MultiQC Edition)
# ------------------------------------------------------------------------------
# Combines:
# - v15.3: Bug fixes, complete 11 signatures, robust dual plasticity testing
# - v14.0: Gene drivers, correlation analysis, rich narrative, BIPLOT
# - v13.0: Multi-FDR TSV exports for MultiQC, simple HTML integration
# - v15.5: ONE unified MultiQC-compatible report with all thresholds
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
# CONFIGURATION
# ==============================================================================
N_TOP_VAR <- opt$n_top_genes
SCORING_METHOD <- opt$scoring
FDR_THRESHOLDS <- c(0.05, 0.01, 0.005, 0.001)
FDR_LABELS <- c("Standard (0.05)", "Strict (0.01)", "Very Strict (0.005)", "Ultra (0.001)")

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================
# Format numbers for display - handles scientific notation for very small p-values
safe_format <- function(x, digits=3) {
    if (is.null(x) || length(x) == 0 || is.na(x)) return("NA")
    if (length(x) > 1) return(paste(sapply(x, safe_format, digits=digits), collapse=", "))
    if (!is.numeric(x)) return(as.character(x))
    if (abs(x) < 0.001 && x != 0) return(formatC(x, format="e", digits=2))
    return(round(x, digits))
}

# Interpret p-value significance levels using standard notation
interpret_p <- function(p) {
    if (is.null(p) || length(p) == 0 || is.na(p)) return("NA")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    if (p < 0.10) return(".")
    return("ns")
}

# Publication-quality theme for plots
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

# Statistical Logger - keeps track of all statistical decisions
stat_log <- list()
add_stat_log <- function(title, content) {
    stat_log[[length(stat_log) + 1]] <<- list(title=title, content=content)
}

# ==============================================================================
# 1. DATA LOADING
# ==============================================================================
cat("LOG [1/10]: Loading Data...\n")

# Create output directory if it doesn't exist
dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)

# Load VST-normalized expression matrix
mat_vst <- as.matrix(read.table(opt$input, header=TRUE, row.names=1, check.names=FALSE))
# Load metadata with sample classifications
meta <- read.csv(opt$meta, row.names=1)

# Align samples between expression matrix and metadata
common <- intersect(colnames(mat_vst), rownames(meta))
if(length(common) < 3) stop("Error: <3 matching samples between VST and Metadata.")
mat_vst <- mat_vst[, common]
meta <- meta[common, , drop=FALSE]

# Ensure Classification column exists
if(!"Classification" %in% colnames(meta)) {
    stop("Metadata must contain a 'Classification' column.")
}

# Set up factor levels for ordered progression (Culture -> Primary -> Recurrent)
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

# Robust color assignment - uses default U251 colors if available, otherwise generates palette
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

# Smart gene ID mapping - detects Ensembl IDs and maps to gene symbols
sample_id <- rownames(mat_vst)[1]
if (grepl("^ENSG", sample_id)) {
    cat("  → Detected Ensembl IDs. Mapping to gene symbols...\n")
    clean_ids <- sub("\\..*", "", rownames(mat_vst))  # Remove version numbers
    symbols <- mapIds(EnsDb.Hsapiens.v86, keys=clean_ids, column="SYMBOL",
                     keytype="GENEID", multiVals="first")

    # Aggregate duplicate gene symbols by taking mean expression
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
cat("LOG [2/10]: Global Structure Analysis...\n")

# Select top variable genes for PCA (reduces noise, focuses on informative features)
top_var <- head(order(apply(mat_sym, 1, var), decreasing=TRUE), N_TOP_VAR)
mat_sig <- mat_sym[top_var, ]

add_stat_log("Feature Selection", sprintf(
    "Selected top %d most variable genes (variance range: %.2f - %.2f)",
    N_TOP_VAR, min(apply(mat_sig, 1, var)), max(apply(mat_sig, 1, var))
))

# PERMANOVA - tests if groups have distinct multivariate expression profiles
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

# PCA with proper dimension handling (fixes bug when N < 10)
pca <- prcomp(t(mat_sig))
n_samples <- ncol(mat_sig)
max_pcs <- min(10, n_samples - 1, ncol(pca$x))  # Can't have more PCs than samples-1

pca_summary <- summary(pca)
n_pcs_to_show <- min(max_pcs, ncol(pca_summary$importance))

var_pc <- round(pca_summary$importance[2, 1:n_pcs_to_show]*100, 1)
cum_var <- round(pca_summary$importance[3, 1:n_pcs_to_show]*100, 1)
pcaData <- data.frame(pca$x[,1:2], Sample=rownames(meta), Class=meta$Classification)

add_stat_log("PCA Variance", sprintf(
    "PC1 = %.1f%%, PC2 = %.1f%%, PC1+PC2 = %.1f%%\n  First %d PCs explain %.1f%% total variance",
    var_pc[1], var_pc[2], cum_var[2], min(5, n_pcs_to_show), cum_var[min(5, n_pcs_to_show)]
))

# PCA Gene Drivers - identify genes with highest loadings on each PC
loadings <- pca$rotation
pc_drivers_text <- ""
pc_drivers_list <- list()

for(i in 1:min(5, ncol(loadings))) {
    pc <- paste0("PC", i)
    # Get top 8 genes with highest absolute loadings
    top <- names(sort(abs(loadings[, i]), decreasing=TRUE)[1:8])
    pc_drivers_list[[pc]] <- top
    pc_drivers_text <- paste0(pc_drivers_text, "\n", pc, ": ", paste(top, collapse=", "))
}

add_stat_log("PCA Gene Drivers", sprintf("Top genes driving each PC:%s", pc_drivers_text))

# Scree Plot - shows variance explained by each PC
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

# PCA Trajectory Plot WITH BIPLOT (gene drivers shown as red arrows)
# Calculate centroids for each group to draw evolutionary trajectory
centroids <- aggregate(cbind(PC1, PC2) ~ Class, data=pcaData, FUN=mean)
centroids <- centroids[match(levels(meta$Classification), centroids$Class), ]
centroids <- na.omit(centroids)

# Create arrows connecting centroids in order (Culture -> Primary -> Recurrent)
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

# BIPLOT: Top gene loadings to show which genes drive PC separation
# Select top 10 genes with highest combined PC1+PC2 loadings
top_genes_load <- loadings[order(sqrt(loadings[,"PC1"]^2 + loadings[,"PC2"]^2), decreasing=TRUE)[1:10], c("PC1", "PC2")]
# Scale gene arrows to fit nicely on plot without overwhelming sample points
gene_arrow_scale <- max(abs(pcaData$PC1)) / max(abs(top_genes_load[,"PC1"])) * 0.8
gene_arrows <- as.data.frame(top_genes_load * gene_arrow_scale)
gene_arrows$Gene <- rownames(gene_arrows)

p_pca <- ggplot(pcaData, aes(x=PC1, y=PC2)) +
    # Trajectory arrows showing evolutionary progression
    {if(nrow(arrow_data) > 0)
        geom_segment(data=arrow_data, aes(x=x, y=y, xend=xend, yend=yend),
                    arrow=arrow(length=unit(0.4,"cm"), type="closed"),
                    color="grey50", linewidth=1.2, inherit.aes=FALSE)
    } +
    # BIPLOT: Gene loading arrows (red) showing gene drivers
    geom_segment(data=gene_arrows, aes(x=0, y=0, xend=PC1, yend=PC2),
                arrow=arrow(length=unit(0.2,"cm")), color="red", alpha=0.4, inherit.aes=FALSE) +
    geom_text_repel(data=gene_arrows, aes(x=PC1, y=PC2, label=Gene),
                   color="red", size=3, segment.alpha=0.4) +
    # Sample points
    geom_point(aes(fill=Class, shape=Class), size=6, color="black", stroke=0.8) +
    geom_text_repel(aes(label=Sample), size=3.5, box.padding=0.5) +
    scale_fill_manual(values=GROUP_COLORS) +
    scale_shape_manual(values=GROUP_SHAPES) +
    labs(title="Evolutionary Trajectory with Gene Drivers",
         subtitle=sprintf("PERMANOVA: R²=%.3f, P=%s | Red arrows = gene loadings", perm_r2, safe_format(perm_p)),
         x=paste0("PC1 (", var_pc[1], "%)"),
         y=paste0("PC2 (", var_pc[2], "%)")) +
    theme_publication()

# Combine PCA and Scree into one composite figure
p_composite <- (p_pca | p_scree) + plot_layout(widths = c(2.5, 1))
ggsave(paste0(opt$out, "_Global_Structure_mqc.png"), p_composite, width=13, height=6)

# Phylogenetic Tree - hierarchical clustering showing sample relationships
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
cat("LOG [3/10]: Scoring Signatures...\n")

sigs <- list()

# Load signatures - either from external file or use built-in comprehensive set
if (!is.null(opt$signatures)) {
    cat("  → Parsing external signature file:", opt$signatures, "\n")
    
    # GMT format parser
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
    # CSV/TSV format parser
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
    # COMPLETE built-in signature set covering 3 major GBM classification systems
    cat("  → Using COMPLETE built-in GBM signatures (Verhaak + Neftel + Garofano)\n")
    sigs <- list(
        # VERHAAK (2010) - TCGA bulk RNA-seq subtypes
        "Verhaak_Classical" = c("EGFR", "NES", "NOTCH3", "JAG1", "HES5", "AKT2", "FGFR2", 
                                "CCND2", "CDK4", "RB1", "SOX2", "SFRP2"),
        "Verhaak_Mesenchymal" = c("CHI3L1", "CD44", "VIM", "RELB", "STAT3", "TRADD", "CASP1",
                                  "TNFRSF1A", "IKBKB", "NFKB1", "MET", "FOSL2", "TIMP1",
                                  "S100A4", "FN1", "COL1A2", "COL3A1", "SERPINE1"),
        "Verhaak_Proneural" = c("OLIG2", "DLL3", "ASCL1", "TCF12", "DCX", "SOX11", "PDGFRA",
                                "NKX2-2", "ERBB3", "NCAM1", "NRCAM", "NKX6-1", "PAX3", "NEUROG1"),
        "Verhaak_Neural" = c("NEFL", "GABRA1", "SYT1", "SLC12A5", "GRIA2", "SYN1", "GABBR1",
                             "SNAP25", "ATP1A3", "STX1A", "GABRG2", "NRXN1", "NLGN1"),
        
        # NEFTEL (2019) - Single-cell GBM states
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
        
        # GAROFANO (2021) - Metabolic subtypes
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

# Calculate gene coverage - what % of signature genes are in the dataset
gene_coverage <- sapply(sigs, function(sig_genes) {
    sum(sig_genes %in% rownames(mat_sym)) / length(sig_genes) * 100
})

add_stat_log("Signature Gene Coverage", sprintf(
    "Mean coverage: %.1f%%\n  %s",
    mean(gene_coverage),
    paste(names(gene_coverage), "=", round(gene_coverage, 1), "%", collapse="\n  ")
))

# Method 1: GSVA - Gene Set Variation Analysis (Hänzelmann et al. 2013)
# Estimates enrichment of gene sets across samples using Gaussian kernel
gsva_res <- suppressWarnings(gsva(mat_sym, sigs, method="gsva", kcdf="Gaussian", verbose=FALSE))

# Method 2: Z-Score Averaging (numerically stable implementation)
# Convert expression to z-scores, then average across signature genes
mat_z <- t(scale(t(mat_sym)))  # Z-score normalization per gene
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

# Compare methods and choose final scoring approach
if(SCORING_METHOD == "both") {
    # Calculate correlation between GSVA and Z-score methods
    cor_methods <- cor(as.vector(gsva_res), as.vector(z_res), use="complete.obs")
    
    add_stat_log("Scoring Method Comparison", sprintf(
        "GSVA vs Z-Score correlation: r = %.3f\n  Interpretation: %s\n  Decision: Using Z-Score (more robust for small gene sets and small N)",
        cor_methods,
        if(cor_methods > 0.8) "High agreement" else if(cor_methods > 0.6) "Moderate agreement" else "Low agreement - interpret with caution"
    ))
    
    # Create agreement plot
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
# 4. PLASTICITY (SHANNON ENTROPY WITH DUAL TESTING)
# ==============================================================================
cat("LOG [4/10]: Plasticity Analysis...\n")

# Calculate Shannon entropy from signature scores
# Higher entropy = more plastic/undifferentiated state (multiple signatures active)
# Lower entropy = more committed state (one dominant signature)
calc_entropy <- function(s) {
    # Clamp extreme values to prevent numerical overflow
    s <- pmin(pmax(s, -10), 10)
    # Convert to probabilities via softmax
    p <- exp(s)/sum(exp(s))
    # Calculate Shannon entropy with small epsilon to prevent log(0)
    -sum(p*log(p + 1e-10), na.rm=TRUE)
}

meta$Plasticity <- apply(final_scores, 2, calc_entropy)

# Smart test selection based on sample size and distribution
n_total <- nrow(meta)
use_nonparametric <- (n_total < 30)  # For small N, Kruskal-Wallis is more powerful

if(use_nonparametric) {
    cat(sprintf("  → Using Kruskal-Wallis (N=%d < 30, more powerful for small samples)\n", n_total))
    plast_test <- kruskal.test(Plasticity ~ Classification, data=meta)
    plast_p <- plast_test$p.value
    plast_stat <- plast_test$statistic
    plast_method <- "Kruskal-Wallis"
    
    add_stat_log("Plasticity: Test Selection", sprintf(
        "Sample size: %d\n  Decision: Kruskal-Wallis (non-parametric)\n  Rationale: More powerful than ANOVA for N<30",
        n_total
    ))
    
    add_stat_log("Plasticity: Group Comparison", sprintf(
        "Kruskal-Wallis H = %.3f, df = %d\n  P = %s %s\n  Interpretation: %s",
        plast_stat, plast_test$parameter,
        safe_format(plast_p), interpret_p(plast_p),
        if(plast_p < 0.05) "Significant differences in plasticity across groups" else "No significant differences"
    ))
} else {
    # For larger samples, test normality first
    shapiro_p <- shapiro.test(meta$Plasticity)$p.value
    
    add_stat_log("Plasticity: Normality Test", sprintf(
        "Shapiro-Wilk W = %.3f, P = %s\n  Interpretation: Data is %s",
        shapiro.test(meta$Plasticity)$statistic,
        safe_format(shapiro_p),
        if(shapiro_p >= 0.05) "normally distributed" else "non-normal"
    ))
    
    if(shapiro_p < 0.05) {
        # Use non-parametric if data is non-normal
        plast_test <- kruskal.test(Plasticity ~ Classification, data=meta)
        plast_p <- plast_test$p.value
        plast_method <- "Kruskal-Wallis"
    } else {
        # Use parametric ANOVA if data is normal
        plast_aov <- aov(Plasticity ~ Classification, data=meta)
        plast_summary <- summary(plast_aov)[[1]]
        plast_p <- plast_summary[["Pr(>F)"]][1]
        plast_f <- plast_summary[["F value"]][1]
        plast_method <- "ANOVA"
        
        add_stat_log("Plasticity: Group Comparison", sprintf(
            "One-way ANOVA F(%d,%d) = %.3f\n  P = %s %s",
            plast_summary[["Df"]][1], plast_summary[["Df"]][2],
            plast_f, safe_format(plast_p), interpret_p(plast_p)
        ))
    }
}

# Calculate descriptive statistics by group
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

# Create plasticity visualization
p_plast <- ggplot(meta, aes(x=Classification, y=Plasticity, fill=Classification)) +
    geom_violin(alpha=0.3, trim=FALSE) +  # Show distribution shape
    geom_boxplot(width=0.2, fill="white", outlier.shape=NA) +  # Show median and quartiles
    geom_jitter(width=0.1, size=3, alpha=0.7) +  # Show individual points
    stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="white") +  # Show mean
    scale_fill_manual(values=GROUP_COLORS) +
    labs(title="Cellular Plasticity (Shannon Entropy)",
         subtitle=sprintf("%s: P=%s %s", plast_method, safe_format(plast_p), interpret_p(plast_p)),
         caption="Diamond = mean, box = median ± IQR, dots = individual samples") +
    theme_publication() +
    theme(legend.position="none")

ggsave(paste0(opt$out, "_Plasticity_mqc.png"), p_plast, width=7, height=6)

# ==============================================================================
# 5. ROBUST LIMMA WITH ARRAY WEIGHTS
# ==============================================================================
cat("LOG [5/10]: Differential Analysis with Sample Quality Weighting...\n")

# Create design matrix for linear model (one column per group)
design <- model.matrix(~0 + meta$Classification)
colnames(design) <- levels(meta$Classification)

# Calculate sample quality weights using arrayWeights
# This down-weights noisy/outlier samples automatically
aw <- arrayWeights(final_scores, design)
weight_df <- data.frame(
    Sample = names(aw),
    Weight = round(aw, 3),
    Group = as.character(meta[names(aw), "Classification"])
)

# Identify potential outliers (samples with weights < mean - 2*SD)
outlier_threshold <- mean(aw) - 2*sd(aw)
outliers <- weight_df$Sample[weight_df$Weight < outlier_threshold]

add_stat_log("Sample Quality Assessment", sprintf(
    "arrayWeights analysis:\n  Mean weight = %.3f\n  SD = %.3f\n  Range = %.3f - %.3f\n  Potential outliers: %s\n  Interpretation: %s",
    mean(aw), sd(aw), min(aw), max(aw),
    if(length(outliers) > 0) paste(outliers, collapse=", ") else "None detected",
    if(length(outliers) > 0) "Some samples show reduced reliability and will be down-weighted" else "All samples show consistent quality"
))

# CRITICAL: Dynamic contrast generation based on actual factor levels
# This ensures the script works with any group names, not just U251 defaults
levs <- levels(meta$Classification)
contrast_formulas <- c()
contrast_names <- c()

# Generate all pairwise comparisons
if(length(levs) >= 2) {
    # Compare second group to first (e.g., Primary vs Culture)
    contrast_formulas <- c(contrast_formulas, sprintf("%s - %s", levs[2], levs[1]))
    contrast_names <- c(contrast_names, sprintf("%s_vs_%s", 
                                                gsub("[^A-Za-z0-9]", "", levs[2]), 
                                                gsub("[^A-Za-z0-9]", "", levs[1])))
}
if(length(levs) >= 3) {
    # Compare third group to first and second (e.g., Recurrent vs Culture, Recurrent vs Primary)
    contrast_formulas <- c(contrast_formulas,
                          sprintf("%s - %s", levs[3], levs[1]),
                          sprintf("%s - %s", levs[3], levs[2]))
    contrast_names <- c(contrast_names,
                       sprintf("%s_vs_%s", gsub("[^A-Za-z0-9]", "", levs[3]), gsub("[^A-Za-z0-9]", "", levs[1])),
                       sprintf("%s_vs_%s", gsub("[^A-Za-z0-9]", "", levs[3]), gsub("[^A-Za-z0-9]", "", levs[2])))
}

cat("\n📋 Contrasts being tested:\n")
for(i in 1:length(contrast_names)) {
    cat(sprintf("  %d. %s: %s\n", i, contrast_names[i], contrast_formulas[i]))
}

# Create contrast matrix
cont.matrix <- makeContrasts(contrasts=contrast_formulas, levels=design)
colnames(cont.matrix) <- contrast_names

# Fit weighted linear model with empirical Bayes moderation
fit <- lmFit(final_scores, design, weights=aw)  # Fit initial model with sample weights
fit <- contrasts.fit(fit, cont.matrix)  # Apply contrasts
fit <- eBayes(fit)  # Empirical Bayes moderation of standard errors

all_coefs <- colnames(cont.matrix)

add_stat_log("Linear Model", sprintf(
    "Contrasts tested: %s\n  Method: Weighted linear regression + empirical Bayes\n  Weights: arrayWeights (automatic quality assessment)\n  Correction: Benjamini-Hochberg FDR",
    paste(contrast_names, collapse=", ")
))

# ==============================================================================
# 6. MULTI-FDR SENSITIVITY ANALYSIS (CRITICAL - THIS IS WHERE YOUR RESULTS ARE!)
# ==============================================================================
cat("LOG [6/10]: Multi-FDR Sensitivity Analysis (generating separate tables for each threshold)...\n")

# Initialize results tracking
sig_counts <- data.frame(FDR_Threshold=FDR_THRESHOLDS, Label=FDR_LABELS)
for(c in all_coefs) sig_counts[[c]] <- 0

# CRITICAL: Generate SEPARATE TSV file for EACH FDR threshold
# This is what MultiQC will parse and display in separate sections
# This is also where your Neftel_AC significance at different thresholds lives!
for(i in 1:length(FDR_THRESHOLDS)) {
    thresh <- FDR_THRESHOLDS[i]
    lbl <- FDR_LABELS[i]
    
    cat(sprintf("\n  → Processing FDR threshold: %s (%.3f)\n", lbl, thresh))
    
    # Create results table for this FDR threshold
    stats_out <- data.frame(Signature=rownames(final_scores))
    
    # Add results for each contrast
    for(coef in all_coefs) {
        tt <- topTable(fit, coef=coef, number=Inf, adjust.method="BH")
        
        # Add statistics columns
        stats_out[[paste0(coef, "_logFC")]] <- round(tt[stats_out$Signature, "logFC"], 3)
        stats_out[[paste0(coef, "_P")]] <- sapply(stats_out$Signature, function(sig) {
            safe_format(tt[sig, "P.Value"])
        })
        stats_out[[paste0(coef, "_adjP")]] <- sapply(stats_out$Signature, function(sig) {
            safe_format(tt[sig, "adj.P.Val"])
        })
        # Mark significant results at this threshold
        stats_out[[paste0(coef, "_Sig")]] <- ifelse(
            tt[stats_out$Signature, "adj.P.Val"] < thresh, "Yes", "No"
        )
        
        # Count significant signatures for summary
        n_sig <- sum(tt[, "adj.P.Val"] < thresh, na.rm=TRUE)
        sig_counts[i, coef] <- n_sig
        
        cat(sprintf("    • %s: %d/%d signatures significant\n", 
                   coef, n_sig, nrow(final_scores)))
    }
    
    # Export to MultiQC-compatible TSV file
    # Format: fdr_<threshold>_stats_mqc.tsv
    thresh_code <- gsub("[^0-9]", "", as.character(thresh))
    tf <- paste0(dirname(opt$out), "/fdr_", thresh_code, "_stats_mqc.tsv")
    
    # Write MultiQC header
    cat(paste0("# id: 'stats_", thresh_code, 
               "'\n# section_name: 'Differential Signatures: ", lbl, 
               "'\n# plot_type: 'table'\n"), file=tf)
    # Write data
    write.table(stats_out, file=tf, sep="\t", quote=FALSE, row.names=FALSE, append=TRUE)
    
    cat(sprintf("    ✓ Exported to: %s\n", basename(tf)))
}

add_stat_log("FDR Sensitivity Analysis", sprintf(
    "Tested %d FDR thresholds across %d contrasts\n  Significant signatures at each threshold:\n  %s",
    length(FDR_THRESHOLDS), length(all_coefs),
    paste(capture.output(print(sig_counts, row.names=FALSE)), collapse="\n  ")
))

# Create FDR sensitivity visualization
sig_counts_long <- sig_counts %>%
    pivot_longer(cols=all_of(all_coefs), names_to="Contrast", values_to="N_Significant") %>%
    mutate(Contrast = factor(Contrast, levels=all_coefs))

p_sensitivity <- ggplot(sig_counts_long, aes(x=Label, y=N_Significant, fill=Contrast)) +
    geom_bar(stat="identity", position="dodge", alpha=0.8) +
    geom_text(aes(label=N_Significant), position=position_dodge(0.9), vjust=-0.5, size=3.5) +
    scale_fill_manual(values=rep(GROUP_COLORS[1:min(3, length(GROUP_COLORS))], length.out=length(all_coefs))) +
    labs(title="FDR Sensitivity Analysis",
         subtitle=sprintf("%d signatures tested across %d contrasts at 4 FDR thresholds", 
                         nrow(final_scores), length(all_coefs)),
         x="FDR Threshold", y="Number of Significant Signatures") +
    theme_publication(base_size=12) +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    ylim(0, nrow(final_scores) + 1)

ggsave(paste0(opt$out, "_FDR_Sensitivity_mqc.png"), p_sensitivity, width=11, height=7)

# Export FDR sensitivity summary for MultiQC
summary_file <- paste0(dirname(opt$out), "/fdr_sensitivity_summary_mqc.tsv")
cat("# id: 'fdr_sensitivity'\n# section_name: 'FDR Sensitivity Summary'\n# plot_type: 'table'\n", 
    file=summary_file)
write.table(sig_counts, file=summary_file, sep="\t", quote=FALSE, row.names=FALSE, append=TRUE)

# ==============================================================================
# 7. UNIFIED TRAJECTORY PLOT (ALL SIGNATURES)
# ==============================================================================
cat("LOG [7/10]: Trajectory Analysis...\n")

# Prepare trajectory data for all signatures
traj_data_list <- list()
for(sig in rownames(final_scores)) {
    df <- data.frame(
        Signature = sig,
        Sample = colnames(final_scores),
        Score = final_scores[sig, ],
        Class = meta$Classification,
        Stage = as.numeric(meta$Classification)  # Convert factor to numeric for correlation
    )
    traj_data_list[[sig]] <- df
}
traj_data <- do.call(rbind, traj_data_list)

# Calculate group means and standard errors for plotting
traj_summary <- traj_data %>%
    group_by(Signature, Class, Stage) %>%
    summarise(Mean = mean(Score), SE = sd(Score)/sqrt(n()), .groups="drop")

# Sort signatures by correlation with stage (increasing to decreasing)
trend_order <- traj_data %>%
    group_by(Signature) %>%
    summarise(Cor = cor(Stage, Score, method="spearman"), .groups="drop") %>%
    arrange(desc(Cor)) %>%
    pull(Signature)

traj_data$Signature <- factor(traj_data$Signature, levels=trend_order)
traj_summary$Signature <- factor(traj_summary$Signature, levels=trend_order)

# Create faceted trajectory plot
p_traj <- ggplot(traj_data, aes(x=Stage, y=Score)) +
    # Confidence ribbon (mean ± SE)
    geom_ribbon(data=traj_summary, aes(x=Stage, ymin=Mean-SE, ymax=Mean+SE, fill=Signature),
                inherit.aes=FALSE, alpha=0.2) +
    # Mean trajectory line (EXPLICIT geom_line ensures all signatures show lines)
    geom_line(data=traj_summary, aes(x=Stage, y=Mean, color=Signature, group=Signature),
              inherit.aes=FALSE, linewidth=1.2, alpha=0.9) +
    # Individual sample points
    geom_point(aes(fill=Class, shape=Class), size=3, alpha=0.6) +
    facet_wrap(~Signature, scales="free_y", ncol=3) +
    scale_x_continuous(breaks=1:length(levs), labels=levs) +
    scale_fill_manual(values=GROUP_COLORS, name="Stage") +
    scale_color_brewer(palette="Set1", guide="none") +
    scale_shape_manual(values=GROUP_SHAPES, name="Stage") +
    labs(title="Signature Trajectories Across Evolution",
         subtitle="Ordered by correlation with stage (top = increasing, bottom = decreasing) | Lines = group means ± SE",
         x="Stage", y="Z-Score Expression") +
    theme_publication(base_size=11) +
    theme(legend.position = "bottom")

ggsave(paste0(opt$out, "_Unified_Trajectories_mqc.png"), p_traj, width=14, height=12)

# Calculate trajectory statistics (correlation with stage)
traj_models <- traj_data %>%
    group_by(Signature) %>%
    summarise(
        Correlation = cor(Stage, Score, method="spearman"),
        Trend_P = cor.test(Stage, Score, method="spearman")$p.value,
        Direction = ifelse(Correlation > 0, "Increasing", "Decreasing"),
        .groups="drop"
    ) %>%
    arrange(desc(Correlation))

add_stat_log("Trajectory Analysis", sprintf(
    "Signature evolution patterns:\n  %s",
    paste(capture.output(print(traj_models, row.names=FALSE)), collapse="\n  ")
))

# ==============================================================================
# 8. HEATMAPS & CORRELATION
# ==============================================================================
cat("LOG [8/10]: Heatmaps & Co-evolution Analysis...\n")

# Signature scores heatmap with annotations
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

# Signature correlation analysis - identifies co-evolving signatures
sig_cor <- cor(t(final_scores), method="pearson")

# Find strong correlations (|r| > 0.6)
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
    add_stat_log("Signature Co-evolution", "No strong correlations (|r| > 0.6) detected")
}

# Correlation heatmap
png(paste0(opt$out, "_Sig_Correlation_mqc.png"), width=8, height=7, units="in", res=300)
Heatmap(sig_cor, name="Pearson\nCorr",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        column_title = "Signature Co-evolution",
        row_names_gp = gpar(fontsize=10),
        column_names_gp = gpar(fontsize=10))
dev.off()

# ==============================================================================
# 9. EXPORT DATA
# ==============================================================================
cat("LOG [9/10]: Exporting Data Files...\n")

write.csv(t(final_scores), paste0(opt$out, "_Scores.csv"))
write.csv(meta, paste0(opt$out, "_Metadata.csv"))
write.csv(weight_df, paste0(opt$out, "_SampleWeights.csv"), row.names=FALSE)
write.csv(plast_stats, paste0(opt$out, "_Plasticity_Statistics.csv"), row.names=FALSE)
write.csv(traj_models, paste0(opt$out, "_Trajectory_Statistics.csv"), row.names=FALSE)

# ==============================================================================
# 10. MULTIQC-COMPATIBLE HTML REPORT
# ==============================================================================
cat("LOG [10/10]: Generating MultiQC-Compatible HTML Report...\n")

summary_html <- paste0(dirname(opt$out), "/analysis_summary_mqc.html")
sink(summary_html)

group_summary <- paste(names(group_sizes), as.vector(group_sizes), sep=": ", collapse="\n")

cat("
<div style='background-color:#f8f9fa; padding:20px; border:2px solid #667eea; border-radius:8px; margin-bottom:30px;'>
    <h2 style='color:#667eea; margin-top:0;'>🧬 Global Subtype Evolution Analysis (v15.5 Ultimate)</h2>
    
    <h3>Dataset Overview</h3>
    <pre style='background:white; padding:10px; border:1px solid #ddd;'>", group_summary, "</pre>
    <p><strong>Scoring Method:</strong> ", scoring_label, "</p>
    ", if(SCORING_METHOD == "both") paste0("<p><strong>Method Agreement (GSVA vs Z-Score):</strong> r = ", round(cor_methods, 3), "</p>") else "", "
    
    <h3>Global Structure</h3>
    <ul style='line-height:1.8;'>
        <li><strong>PERMANOVA:</strong> F = ", sprintf("%.3f", perm_f), 
            ", R² = ", sprintf("%.3f", perm_r2), " (", sprintf("%.1f", perm_r2*100), "% variance explained)",
            ", P = ", safe_format(perm_p), " ", interpret_p(perm_p), "</li>
        <li><strong>Plasticity (", plast_method, "):</strong> P = ", safe_format(plast_p), " ", interpret_p(plast_p), "</li>
    </ul>
    
    <h3>Sample Quality Weights</h3>
    <p style='color:#666; font-size:13px;'>arrayWeights automatically identifies and down-weights noisy samples. Lower weights indicate reduced reliability.</p>
    <table style='border-collapse: collapse; width:100%; margin:15px 0;'>
        <tr style='background:#667eea; color:white;'>
            <th style='padding:8px; border:1px solid #555;'>Sample</th>
            <th style='padding:8px; border:1px solid #555;'>Weight</th>
            <th style='padding:8px; border:1px solid #555;'>Group</th>
        </tr>
")

for(j in order(weight_df$Weight)) {
    # Highlight potential outliers in red
    row_color <- if(weight_df$Weight[j] < mean(aw) - sd(aw)) "#ffe6e6" else "white"
    cat(sprintf("        <tr style='background:%s;'><td style='padding:8px; border:1px solid #ccc;'>%s</td>
            <td style='padding:8px; border:1px solid #ccc;'>%.3f</td>
            <td style='padding:8px; border:1px solid #ccc;'>%s</td></tr>\n",
            row_color, weight_df$Sample[j], weight_df$Weight[j], weight_df$Group[j]))
}

cat("    </table>
    
    <h3>PCA Gene Drivers (Biplot Arrows)</h3>
    <p style='color:#666; font-size:13px;'>Top genes with highest loadings on each principal component (red arrows in PCA plot).</p>
    <pre style='font-size:11px; background:white; padding:10px; border:1px solid #ddd;'>", pc_drivers_text, "</pre>
    
    <h3>Plasticity Statistics by Group</h3>
    <p style='color:#666; font-size:13px;'>Shannon entropy calculated from signature scores. Higher = more plastic/undifferentiated.</p>
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
    
    <h3>Signature Trajectory Analysis</h3>
    <p style='color:#666; font-size:13px;'>Spearman correlation of signature scores with evolutionary stage (1=Culture, 2=Primary, 3=Recurrent).</p>
    <table style='border-collapse: collapse; width:100%; margin:15px 0;'>
        <tr style='background:#667eea; color:white;'>
            <th style='padding:8px; border:1px solid #555;'>Signature</th>
            <th style='padding:8px; border:1px solid #555;'>Correlation (ρ)</th>
            <th style='padding:8px; border:1px solid #555;'>P-value</th>
            <th style='padding:8px; border:1px solid #555;'>Direction</th>
            <th style='padding:8px; border:1px solid #555;'>Significant</th>
        </tr>
")

for(i in 1:nrow(traj_models)) {
    sig_marker <- if(traj_models$Trend_P[i] < 0.05) {
        "<span style='color:#28a745; font-weight:bold;'>✓ Yes</span>"
    } else {
        "<span style='color:#6c757d;'>— No</span>"
    }
    
    cat(sprintf("        <tr><td style='padding:8px; border:1px solid #ccc;'><strong>%s</strong></td>
            <td style='padding:8px; border:1px solid #ccc;'>%.3f</td>
            <td style='padding:8px; border:1px solid #ccc;'>%s</td>
            <td style='padding:8px; border:1px solid #ccc;'>%s</td>
            <td style='padding:8px; border:1px solid #ccc;'>%s</td></tr>\n",
            traj_models$Signature[i], traj_models$Correlation[i],
            safe_format(traj_models$Trend_P[i]), traj_models$Direction[i], sig_marker))
}

cat("    </table>
    
    <h3>FDR Sensitivity Analysis Results</h3>
    <p style='color:#666; font-size:13px;'><strong>THIS IS WHERE YOUR NEFTEL_AC RESULTS ARE!</strong> 
    Number of significant signatures at each FDR threshold. Detailed tables for each threshold are in separate MultiQC sections above.</p>
    <table style='border-collapse: collapse; width:100%; margin:15px 0;'>
        <tr style='background:#667eea; color:white;'>
            <th style='padding:8px; border:1px solid #555;'>Threshold</th>")

for(coef in all_coefs) {
    cat(sprintf("<th style='padding:8px; border:1px solid #555;'>%s</th>", gsub("_", " ", coef)))
}
cat("</tr>\n")

for(i in 1:nrow(sig_counts)) {
    cat(sprintf("        <tr><td style='padding:8px; border:1px solid #ccc;'><strong>%s</strong></td>", 
                sig_counts$Label[i]))
    for(coef in all_coefs) {
        count_val <- sig_counts[i, coef]
        # Highlight non-zero counts
        bg_color <- if(count_val > 0) "#e7f3ff" else "white"
        cat(sprintf("<td style='padding:8px; border:1px solid #ccc; text-align:center; background:%s;'><strong>%d</strong></td>", 
                   bg_color, count_val))
    }
    cat("</tr>\n")
}

cat("    </table>
    <p style='font-size:12px; color:#666; margin-top:10px; padding:10px; background:#fff3cd; border-left:4px solid #ffc107;'>
    <strong>📌 Important:</strong> Check the separate <strong>\"Differential Signatures\"</strong> sections above in MultiQC for detailed results at each FDR threshold. 
    These tables show logFC, P-values, and adjusted P-values for all signatures at 0.05, 0.01, 0.005, and 0.001 FDR cutoffs.
    </p>
</div>

<div style='background-color:#fff; padding:15px; border:1px solid #ddd; border-radius:5px; margin-top:20px;'>
    <h3>Statistical Methods</h3>
    <ul style='font-size:13px; line-height:1.8;'>
        <li><strong>Scoring:</strong> ", scoring_label, " (", 
            if(SCORING_METHOD == "both") paste0("GSVA vs Z-Score correlation: r = ", round(cor_methods, 3)) else "user-specified", 
            ")</li>
        <li><strong>Model:</strong> limma with arrayWeights + empirical Bayes moderation</li>
        <li><strong>Contrasts:</strong> ", paste(contrast_names, collapse=", "), "</li>
        <li><strong>FDR Correction:</strong> Benjamini-Hochberg (4 thresholds tested: 0.05, 0.01, 0.005, 0.001)</li>
        <li><strong>Plasticity Test:</strong> ", plast_method, " (automatic selection based on N and normality)</li>
        <li><strong>Trajectory Analysis:</strong> Spearman correlation with evolutionary stage</li>
        <li><strong>Software:</strong> R ", as.character(R.version.string), 
        ", limma ", as.character(packageVersion("limma")),
        ", GSVA ", as.character(packageVersion("GSVA")), "</li>
    </ul>
    
    <h3>Interpretation Guide</h3>
    <ul style='font-size:13px; line-height:1.8;'>
        <li><strong>*** (P < 0.001):</strong> Highly significant - very strong evidence</li>
        <li><strong>** (P < 0.01):</strong> Very significant - strong evidence</li>
        <li><strong>* (P < 0.05):</strong> Significant - moderate evidence</li>
        <li><strong>. (P < 0.10):</strong> Trend - suggestive but not conclusive</li>
        <li><strong>ns (P ≥ 0.10):</strong> Not significant - insufficient evidence</li>
    </ul>
    
    <h3>Key Metrics</h3>
    <ul style='font-size:13px; line-height:1.8;'>
        <li><strong>R² (PERMANOVA):</strong> Proportion of variance explained by groups (0-1 scale)</li>
        <li><strong>Shannon Entropy:</strong> Measure of cellular plasticity; higher = more plastic state</li>
        <li><strong>logFC:</strong> Log2 fold change; positive = higher in second group, negative = higher in first</li>
        <li><strong>Sample Weights:</strong> Quality metric from arrayWeights; ~1.0 is normal, <0.5 indicates noise</li>
    </ul>
</div>

<div style='text-align:center; margin-top:30px; padding:15px; background:#f0f0f0; border-radius:5px;'>
    <p style='margin:0; color:#666; font-size:12px;'>
    Generated by <strong>U251 Global Subtype Analysis Pipeline v15.5 (Ultimate MultiQC Edition)</strong><br>
    Analysis timestamp: ", as.character(Sys.time()), "<br>
    Random seed: 12345 | Top variable genes: ", N_TOP_VAR, "
    </p>
</div>
")

sink()

# Session info
writeLines(capture.output(sessionInfo()), paste0(dirname(opt$out), "/sessionInfo.txt"))

# ==============================================================================
# COMPLETION SUMMARY
# ==============================================================================
cat("\n╔═══════════════════════════════════════════════╗\n")
cat("║      ANALYSIS COMPLETE - v15.5 ULTIMATE       ║\n")
cat("╚═══════════════════════════════════════════════╝\n\n")
cat(sprintf("✓ Analyzed %d samples, %d genes, %d signatures\n", ncol(mat_vst), nrow(mat_sym), nrow(final_scores)))
cat(sprintf("✓ Generated %d FDR-specific MultiQC tables\n", length(FDR_THRESHOLDS)))
cat(sprintf("✓ PCA with BIPLOT showing gene drivers\n"))
cat(sprintf("✓ Unified MultiQC-compatible HTML report\n"))
cat(sprintf("\nKey features:\n"))
cat(sprintf("  • Multi-FDR sensitivity (%s)\n", paste(FDR_LABELS, collapse=", ")))
cat(sprintf("  • Complete 11 signatures (Verhaak: 4, Neftel: 4, Garofano: 3)\n"))
cat(sprintf("  • Robust dual plasticity testing (smart test selection)\n"))
cat(sprintf("  • PCA with gene drivers (biplot)\n"))
cat(sprintf("  • Gene co-evolution analysis\n"))
cat(sprintf("  • Dynamic contrast generation\n"))
cat(sprintf("  • SINGLE MultiQC-integrated report\n"))
cat(sprintf("\nMultiQC will find:\n"))
cat(sprintf("  • analysis_summary_mqc.html (main report with all summaries)\n"))
cat(sprintf("  • fdr_005_stats_mqc.tsv (FDR < 0.05 detailed results)\n"))
cat(sprintf("  • fdr_001_stats_mqc.tsv (FDR < 0.01 detailed results)\n"))
cat(sprintf("  • fdr_0005_stats_mqc.tsv (FDR < 0.005 detailed results)\n"))
cat(sprintf("  • fdr_0001_stats_mqc.tsv (FDR < 0.001 detailed results)\n"))
cat(sprintf("  • fdr_sensitivity_summary_mqc.tsv (summary table)\n"))
cat(sprintf("\n📊 FDR Sensitivity Summary:\n"))
print(sig_counts)
cat(sprintf("\n🎯 Your Neftel_AC results should be in the FDR tables above!\n"))
cat(sprintf("   Check each threshold to see where it becomes significant.\n\n"))
