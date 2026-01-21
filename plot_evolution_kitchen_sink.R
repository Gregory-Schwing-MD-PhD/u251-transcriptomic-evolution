#!/usr/bin/env Rscript
#!/usr/bin/env Rscript

# ADDED: Global reproducibility
set.seed(12345)

# ------------------------------------------------------------------------------
# EVOLUTIONARY KITCHEN SINK v16.1 (Publication + AI) - EXPANDED PLOTS
# ------------------------------------------------------------------------------
# - Enhanced error handling and validation
# - Configurable parameters as constants
# - Better logging and diagnostics
# - Additional plots: Heatmaps, MA plots, GSEA running scores, Network plots
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(ape); library(ggrepel)
    library(EnsDb.Hsapiens.v86); library(clusterProfiler); library(enrichplot)
    library(GSVA); library(GSEABase); library(ComplexHeatmap); library(circlize)
    library(EnhancedVolcano); library(data.table); library(tidyr)
    if (requireNamespace("ggupset", quietly = TRUE)) library(ggupset)
})

# ==============================================================================
# CONFIGURATION CONSTANTS
# ==============================================================================
N_TOP_VARIABLE_GENES <- 500
N_TOP_DE_HEATMAP <- 50
N_TOP_BOXPLOT_GENES <- 6
GSEA_SHOW_CATEGORIES_DOT <- 15
GSEA_SHOW_CATEGORIES_EMAP <- 20
GSEA_SHOW_CATEGORIES_RUNNING <- 3
GSEA_SHOW_CATEGORIES_CNET <- 5
PADJ_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 1.0
VOLCANO_POINT_SIZE <- 2.0
VOLCANO_LABEL_SIZE <- 3.0
VOLCANO_MAX_OVERLAPS <- 30
N_TOP_DE_GENES_EXPORT <- 50
N_TOP_GSEA_EXPORT <- 50

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

# Map gene IDs to symbols with consistent handling
map_genes_to_symbols <- function(gene_ids, db = EnsDb.Hsapiens.v86) {
    clean_ids <- sub("\\..*", "", gene_ids)
    symbols <- mapIds(db, keys=clean_ids, column="SYMBOL",
                     keytype="GENEID", multiVals="first")
    ifelse(is.na(symbols), gene_ids, symbols)
}

# Save plots with comprehensive error handling
save_plot <- function(plot_obj, filename_base, w=9, h=8) {
    tryCatch({
        if (inherits(plot_obj, "Heatmap")) {
            pdf(paste0(filename_base, ".pdf"), width=w, height=h)
            draw(plot_obj)
            dev.off()
            png(paste0(filename_base, ".png"), width=w, height=h, units="in", res=300)
            draw(plot_obj)
            dev.off()
        } else if (inherits(plot_obj, "gg")) {
            ggsave(paste0(filename_base, ".pdf"), plot_obj, width=w, height=h)
            ggsave(paste0(filename_base, ".png"), plot_obj, width=w, height=h, dpi=300, bg="white")
        } else {
            # For base R plots and other objects
            pdf(paste0(filename_base, ".pdf"), width=w, height=h)
            print(plot_obj)
            dev.off()
            png(paste0(filename_base, ".png"), width=w, height=h, units="in", res=300)
            print(plot_obj)
            dev.off()
        }
        cat(paste0("SUCCESS: Saved ", basename(filename_base), "\n"))
        return(TRUE)
    }, error = function(e) {
        cat(paste0("ERROR: Failed to save ", basename(filename_base), ": ", e$message, "\n"))
        return(FALSE)
    })
}

# ==============================================================================
# ARGUMENT PARSING & SETUP
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    stop("Usage: script.R <results_dir> <out_prefix> <meta> <contrasts> <gmt>")
}

results_dir   <- args[1]
output_prefix <- args[2]
meta_file     <- args[3]
contrast_file <- args[4]
gmt_file      <- args[5]

# Validate input files
if (!dir.exists(results_dir)) stop("FATAL: Results directory not found: ", results_dir)
if (!file.exists(meta_file)) stop("FATAL: Metadata file not found: ", meta_file)
if (!file.exists(contrast_file)) stop("FATAL: Contrast file not found: ", contrast_file)
if (!file.exists(gmt_file)) stop("FATAL: GMT file not found: ", gmt_file)

# Create output directories
dir.create(dirname(output_prefix), showWarnings = FALSE, recursive = TRUE)
gemini_dir <- paste0(dirname(output_prefix), "/GEMINI_DATA")
dir.create(gemini_dir, showWarnings = FALSE)

# ==============================================================================
# 1. LOAD AND VALIDATE DATA
# ==============================================================================
cat("LOG [1/6]: Loading and Validating Data...\n")

vst_file <- file.path(results_dir, "tables/processed_abundance/all.vst.tsv")
if (!file.exists(vst_file)) stop("FATAL: Pipeline VST file not found: ", vst_file)

mat_vst <- tryCatch({
    as.matrix(read.table(vst_file, header=TRUE, row.names=1, check.names=FALSE))
}, error = function(e) {
    stop("FATAL: Failed to read VST file: ", e$message)
})

# Validate VST matrix
if (nrow(mat_vst) == 0) stop("FATAL: VST matrix is empty")
if (ncol(mat_vst) < 3) stop("FATAL: Need at least 3 samples, found ", ncol(mat_vst))
cat(paste0("  Matrix dimensions: ", nrow(mat_vst), " genes x ", ncol(mat_vst), " samples\n"))
cat(paste0("  Memory usage: ", format(object.size(mat_vst), units="Mb"), "\n"))

# Load and validate metadata
meta <- tryCatch({
    read.csv(meta_file, row.names=1)
}, error = function(e) {
    stop("FATAL: Failed to read metadata file: ", e$message)
})

# Match samples
common <- intersect(colnames(mat_vst), rownames(meta))
if (length(common) == 0) stop("FATAL: No sample overlap between VST matrix and metadata")
if (length(common) < ncol(mat_vst)) {
    cat(paste0("  WARNING: ", ncol(mat_vst) - length(common), " samples in VST not found in metadata\n"))
}

mat_vst <- mat_vst[, common]
# ADDED: Defensive sample order assertion
stopifnot(identical(colnames(mat_vst), rownames(meta)))

meta <- meta[common, , drop=FALSE]
cat(paste0("  Matched ", length(common), " samples\n"))

# Set Classification factor levels if present
if("Classification" %in% colnames(meta)) {
    meta$Classification <- factor(meta$Classification,
                                  levels=c("Culture_U2", "Primary_U2", "Recurrent_U2"))
}

# Load GMT file
gmt <- tryCatch({
    read.gmt(gmt_file)
}, error = function(e) {
    stop("FATAL: Failed to parse GMT file: ", e$message)
})
if (nrow(gmt) == 0) stop("FATAL: GMT file contains no gene sets")
cat(paste0("  Loaded ", length(unique(gmt$term)), " gene sets from GMT\n"))

# ==============================================================================
# 2. EVOLUTIONARY TREE
# ==============================================================================
cat("LOG [2/6]: Generating Phylogenetic Tree...\n")

top_var <- head(order(apply(mat_vst, 1, var), decreasing=TRUE), N_TOP_VARIABLE_GENES)
mat_sig <- mat_vst[top_var, ]

tryCatch({
    phylo_tree <- as.phylo(hclust(dist(t(mat_sig)), method="ward.D2"))

    if("C2B" %in% phylo_tree$tip.label) {
        rooted_tree <- root(phylo_tree, outgroup="C2B", resolve.root=TRUE)
    } else {
        rooted_tree <- phylo_tree
    }

    pdf(paste0(output_prefix, "_Phylogenetic_Tree_mqc.pdf"), width=9, height=8)
    plot(rooted_tree, main="", type="phylogram", edge.width=2, cex=1.0)
    dev.off()

    cat("  SUCCESS: Phylogenetic tree saved\n")
}, error = function(e) {
    cat(paste0("  ERROR: Tree generation failed: ", e$message, "\n"))
})

# ==============================================================================
# 3. PCA WITH GBM SUBTYPE CLASSIFICATION
# ==============================================================================
cat("LOG [3/6]: Generating PCA with Subtype Classification...\n")

# Map gene IDs to symbols using utility function
mapped_syms <- map_genes_to_symbols(rownames(mat_vst))
# ADDED: Log duplicated gene symbols
dup_n <- sum(duplicated(mapped_syms))
cat(paste0("  Collapsed ", dup_n, " duplicated gene symbols by averaging\n"))

mat_vst_sym <- mat_vst
rownames(mat_vst_sym) <- mapped_syms

# Handle duplicates by averaging
mat_vst_sym <- as.data.frame(mat_vst_sym) %>%
    tibble::rownames_to_column("symbol") %>%
    filter(!is.na(symbol)) %>%
    group_by(symbol) %>%
    summarise(across(everything(), mean)) %>%
    tibble::column_to_rownames("symbol") %>%
    as.matrix()

cat(paste0("  Mapped to ", nrow(mat_vst_sym), " unique gene symbols\n"))

# GBM molecular signatures
gbm_sigs <- list(
    Mesenchymal = c("CHI3L1", "CD44", "VIM", "RELB", "STAT3", "MET", "TRADD", "MMP9", "TIMP1"),
    Classical   = c("EGFR", "AKT2", "NOTCH3", "JAG1", "CCND2", "F3", "PDGFA", "NES"),
    Proneural   = c("PDGFRA", "IDH1", "OLIG2", "SOX2", "NKX2-2", "OLIG1", "TP53"),
    Neural      = c("NEFL", "GABRA1", "SYT1", "SLC12A5", "MBP", "GABRG2")
)

# GSVA subtype scores
gsva_sub <- tryCatch({
    gsva(mat_vst_sym, gbm_sigs, method="gsva", kcdf="Gaussian", verbose=FALSE)
}, error = function(e) {
    stop("FATAL: GSVA failed: ", e$message)
})

gsva_norm <- apply(gsva_sub, 2, function(x) {
    round((x - min(x) + 0.1) / sum(x - min(x) + 0.1) * 100, 0)
})

# PCA calculation
pca <- prcomp(t(mat_sig))
pca_var <- round(summary(pca)$importance[2, 1:2] * 100, 1)

pcaData <- data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    Sample = rownames(meta),
    Class = meta$Classification
)
pcaData$Dominant_Subtype <- apply(gsva_sub, 2, function(x) rownames(gsva_sub)[which.max(x)])

# Export PCA coordinates
write.csv(pcaData, paste0(gemini_dir, "/pca_coordinates.csv"), row.names=FALSE)

# Calculate centroids
centroids <- aggregate(cbind(PC1, PC2) ~ Class, data=pcaData, FUN=mean)
centroids <- centroids[order(centroids$Class),]

# HARDCODED label positions (as requested)
centroids$Label_X <- NA
centroids$Label_Y <- NA
centroids[centroids$Class == "Culture_U2", c("Label_X", "Label_Y")] <- c(-20, 5)
centroids[centroids$Class == "Primary_U2", c("Label_X", "Label_Y")] <- c(-5, 5)
centroids[centroids$Class == "Recurrent_U2", c("Label_X", "Label_Y")] <- c(-5, -10)

# Compute subtype composition by class
comp_summ <- as.data.frame(t(gsva_norm))
comp_summ$Class <- meta[rownames(comp_summ), "Classification"]
comp_summ <- comp_summ %>%
    group_by(Class) %>%
    summarise(across(everything(), mean))

# Create centroid labels
centroids$Label <- apply(centroids, 1, function(row) {
    cls <- row["Class"]
    vals <- comp_summ[comp_summ$Class == cls, -1]
    pcts <- paste(names(vals), paste0(round(as.numeric(vals),0), "%"), sep=": ", collapse="\n")
    paste0(sub("_U2", "", cls), "\n", pcts)
})

# Arrow data for evolutionary trajectory
arrow_data <- data.frame(
    x_start = centroids$PC1[-nrow(centroids)],
    y_start = centroids$PC2[-nrow(centroids)],
    x_end = centroids$PC1[-1],
    y_end = centroids$PC2[-1]
)

# Subtype colors
subtype_colors <- c(
    "Mesenchymal" = "#E41A1C",
    "Proneural" = "#377EB8",
    "Classical" = "#4DAF4A",
    "Neural" = "#984EA3",
    "Other/Unclassified" = "grey70"
)

# Create PCA plot
p_pca <- ggplot(pcaData, aes(x=PC1, y=PC2)) +
    geom_segment(data=arrow_data,
                 aes(x=x_start, y=y_start, xend=x_end, yend=y_end),
                 arrow=arrow(length=unit(0.4,"cm"), type="closed"),
                 color="grey50", linewidth=1.5, inherit.aes=FALSE) +
    geom_point(aes(fill=Dominant_Subtype, shape=Class),
               size=7, color="black", stroke=0.8) +
    geom_label(data=centroids, aes(x=Label_X, y=Label_Y, label=Label),
               fill="white", alpha=0.9, size=3.5, fontface="plain",
               inherit.aes=FALSE) +
    scale_fill_manual(values=subtype_colors) +
    scale_shape_manual(values=c(21, 24, 22),
                       labels=function(x) sub("_U2", "", x)) +
    guides(
        shape = guide_legend(override.aes = list(fill = "black")),
        fill = guide_legend(override.aes = list(shape = 21, size = 5))
    ) +
    theme_bw(base_size=14) +
    theme(
        legend.position="right",
        legend.spacing.y = unit(0.8, 'cm'),
        legend.key.height = unit(1.2, 'cm'),
        panel.grid.minor = element_blank()
    ) +
    labs(
        x = paste0("PC1 (", pca_var[1], "%)"),
        y = paste0("PC2 (", pca_var[2], "%)"),
        title = NULL,
        subtitle = NULL
    )

save_plot(p_pca, paste0(output_prefix, "_PCA_Publication_mqc"), 9, 8)

# ==============================================================================
# 4. SAMPLE QC HEATMAPS
# ==============================================================================
cat("LOG [4/6]: Generating QC Heatmaps...\n")

# Sample-to-sample correlation heatmap
tryCatch({
    sample_cor <- cor(mat_vst, method="pearson")

    ht_cor <- Heatmap(
        sample_cor,
        name = "Pearson\nCorrelation",
        column_split = meta$Classification,
        row_split = meta$Classification,
        col = colorRamp2(c(0.7, 0.85, 1), c("blue", "white", "red")),
        show_row_names = FALSE,
        show_column_names = FALSE,
        border = TRUE,
        column_title = NULL,
        row_title = NULL
    )

    save_plot(ht_cor, paste0(output_prefix, "_SampleCorrelation_mqc"), 8, 8)
}, error = function(e) {
    cat(paste0("  ERROR: Sample correlation heatmap failed: ", e$message, "\n"))
})

# Subtype signature heatmap
tryCatch({
    ht_subtype <- Heatmap(
        gsva_norm,
        name = "Subtype\nScore (%)",
        column_split = meta$Classification,
        cluster_columns = TRUE,
        cluster_rows = FALSE,
        show_column_names = FALSE,
        col = colorRamp2(c(0, 50, 100), c("blue", "white", "red")),
        border = TRUE,
        column_title = NULL,
        row_title = NULL
    )

    save_plot(ht_subtype, paste0(output_prefix, "_SubtypeSignatures_mqc"), 9, 5)
}, error = function(e) {
    cat(paste0("  ERROR: Subtype signature heatmap failed: ", e$message, "\n"))
})

# ==============================================================================
# 5. CONTRAST-SPECIFIC VISUALIZATIONS
# ==============================================================================
cat("LOG [5/6]: Processing Differential Expression Contrasts...\n")

contrasts <- tryCatch({
    read.csv(contrast_file)
}, error = function(e) {
    stop("FATAL: Failed to read contrasts file: ", e$message)
})

if (!"id" %in% colnames(contrasts)) {
    stop("FATAL: Contrasts file must contain 'id' column")
}

prompt_summary <- list()
processed_count <- 0
failed_count <- 0

for(i in 1:nrow(contrasts)) {
    cid <- contrasts$id[i]
    cat(paste0("  Processing contrast: ", cid, "\n"))

    res_file <- file.path(results_dir, "tables/differential",
                         paste0(cid, ".deseq2.results.tsv"))

    if(!file.exists(res_file)) {
        cat(paste0("    WARNING: Results file not found, skipping\n"))
        failed_count <- failed_count + 1
        next
    }

    tryCatch({
        # Load DESeq2 results
        res_df <- read.table(res_file, header=TRUE, sep="\t", quote="")

        # Fix rownames if numeric
        if (grepl("^[0-9]+$", rownames(res_df)[1])) {
            id_col_idx <- which(colnames(res_df) %in%
                              c("gene_id", "gene_name", "id", "GeneID"))
            if(length(id_col_idx) > 0) {
                rownames(res_df) <- res_df[, id_col_idx[1]]
            } else {
                rownames(res_df) <- res_df[,1]
            }
        }

        # Add gene symbols if missing
        if(!"symbol" %in% colnames(res_df)) {
            res_df$symbol <- map_genes_to_symbols(rownames(res_df))
        }

        # Validate essential columns
        if(!all(c("log2FoldChange", "padj") %in% colnames(res_df))) {
            stop("Missing required columns: log2FoldChange or padj")
        }

        # Export top DE genes
        top_de <- res_df %>%
            filter(!is.na(padj)) %>%
            arrange(padj) %>%
            head(N_TOP_DE_GENES_EXPORT)
        write.csv(top_de, paste0(gemini_dir, "/DE_Genes_", cid, ".csv"),
                 row.names=TRUE)

        # ---- VOLCANO PLOT ----
        p_vol <- EnhancedVolcano(
            res_df,
            lab = res_df$symbol,
            x = 'log2FoldChange',
            y = 'padj',
            title = NULL,
            subtitle = NULL,
            caption = NULL,
            pCutoff = PADJ_CUTOFF,
            FCcutoff = LOG2FC_CUTOFF,
            pointSize = VOLCANO_POINT_SIZE,
            labSize = VOLCANO_LABEL_SIZE,
            drawConnectors = TRUE,
            max.overlaps = VOLCANO_MAX_OVERLAPS
        )
        save_plot(p_vol, paste0(output_prefix, "_Volcano_", cid, "_mqc"), 9, 8)

        # ---- MA PLOT ----
        tryCatch({
            if("baseMean" %in% colnames(res_df)) {
                ma_data <- res_df %>%
                    filter(!is.na(padj), !is.na(baseMean), is.finite(log2FoldChange)) %>%
                    mutate(sig = padj < PADJ_CUTOFF & abs(log2FoldChange) > LOG2FC_CUTOFF)

                p_ma <- ggplot(ma_data, aes(x=log10(baseMean + 1), y=log2FoldChange, color=sig)) +
                    geom_point(alpha=0.4, size=1) +
                    scale_color_manual(values=c("grey60", "red"),
                                     labels=c("Not Sig", "Significant")) +
                    geom_hline(yintercept=c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF),
                             linetype="dashed", color="blue", alpha=0.6) +
                    geom_hline(yintercept=0, color="black") +
                    theme_bw(base_size=14) +
                    theme(legend.position="top") +
                    labs(
                        x = "log10(Mean Expression)",
                        y = "log2(Fold Change)",
                        color = NULL,
                        title = NULL
                    )

                save_plot(p_ma, paste0(output_prefix, "_MA_", cid, "_mqc"), 9, 8)
            }
        }, error = function(e) {
            cat(paste0("    WARNING: MA plot failed: ", e$message, "\n"))
        })

        # ---- HEATMAP OF TOP DE GENES ----
        tryCatch({
            top_genes <- res_df %>%
                filter(!is.na(padj)) %>%
                arrange(padj) %>%
                head(N_TOP_DE_HEATMAP) %>%
                pull(symbol)

            # Filter to genes present in mat_vst_sym
            top_genes <- intersect(top_genes, rownames(mat_vst_sym))

            if(length(top_genes) >= 10) {
                heatmap_mat <- mat_vst_sym[top_genes, ]
                heatmap_mat <- t(scale(t(heatmap_mat)))  # Z-score normalization

                ht_de <- Heatmap(
                    heatmap_mat,
                    name = "Z-score",
                    column_split = meta$Classification,
                    show_column_names = FALSE,
                    cluster_columns = TRUE,
                    cluster_rows = TRUE,
                    show_row_names = TRUE,
                    row_names_gp = gpar(fontsize = 8),
                    border = TRUE,
                    column_title = NULL,
                    row_title = NULL
                )

                save_plot(ht_de, paste0(output_prefix, "_Heatmap_TopDE_", cid, "_mqc"), 9, 10)
            }
        }, error = function(e) {
            cat(paste0("    WARNING: DE heatmap failed: ", e$message, "\n"))
        })

        # ---- BOXPLOTS OF TOP GENES ----
        tryCatch({
            top_box_genes <- res_df %>%
                filter(!is.na(padj)) %>%
                arrange(padj) %>%
                head(N_TOP_BOXPLOT_GENES) %>%
                pull(symbol)

            top_box_genes <- intersect(top_box_genes, rownames(mat_vst_sym))

            if(length(top_box_genes) >= 3) {
                box_data <- data.frame(
                    expression = as.vector(mat_vst_sym[top_box_genes, ]),
                    gene = rep(top_box_genes, each = ncol(mat_vst_sym)),
                    classification = rep(meta$Classification, times = length(top_box_genes))
                )

                p_box <- ggplot(box_data, aes(x = classification, y = expression, fill = classification)) +
                    geom_boxplot(outlier.size = 1) +
                    facet_wrap(~gene, scales = "free_y", ncol = 3) +
                    theme_bw(base_size=12) +
                    theme(
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        legend.position = "none",
                        strip.text = element_text(face = "bold")
                    ) +
                    labs(
                        x = NULL,
                        y = "VST Expression",
                        title = NULL
                    ) +
                    scale_fill_brewer(palette = "Set2")

                save_plot(p_box, paste0(output_prefix, "_Boxplot_TopGenes_", cid, "_mqc"), 10, 8)
            }
        }, error = function(e) {
            cat(paste0("    WARNING: Boxplot failed: ", e$message, "\n"))
        })

        # ---- GSEA ----
# ADDED: GSEA ranking metric documentation
# NOTE: GSEA ranking metric = log2FoldChange (directional biology-first choice)

        # Prepare gene list with robust filtering
        gene_list <- res_df %>%
            filter(!is.na(log2FoldChange), !is.na(symbol),
                   is.finite(log2FoldChange)) %>%
            distinct(symbol, .keep_all = TRUE) %>%
            arrange(desc(log2FoldChange)) %>%
            pull(log2FoldChange, name = symbol)

        if (length(gene_list) < 10) {
            cat("    WARNING: Too few genes for GSEA, skipping\n")
            next
        }

        gsea_out <- GSEA(
            gene_list,
            TERM2GENE = gmt,
            pvalueCutoff = 1,
            verbose = FALSE,
            eps = 1e-50
        )

        if(!is.null(gsea_out) && nrow(gsea_out) > 0) {
            # Export GSEA results
            gsea_export <- gsea_out@result %>%
                arrange(p.adjust) %>%
                head(N_TOP_GSEA_EXPORT) %>%
                select(ID, NES, p.adjust, core_enrichment)
            write.csv(gsea_export,
                     paste0(gemini_dir, "/GSEA_Results_", cid, ".csv"),
                     row.names=FALSE)

            # Generate prompt summary
            top_up <- head(gsea_out@result[gsea_out@result$NES > 0, "ID"], 5)
            top_dn <- head(gsea_out@result[gsea_out@result$NES < 0, "ID"], 5)
            prompt_summary[[cid]] <- paste0(
                "Contrast ", cid, ":\n",
                "- Top Up: ", paste(top_up, collapse=", "), "\n",
                "- Top Dn: ", paste(top_dn, collapse=", "), "\n"
            )

            # Calculate term similarity
            gsea_out <- pairwise_termsim(gsea_out)

            # Dotplot
            p_dot <- dotplot(gsea_out, showCategory=GSEA_SHOW_CATEGORIES_DOT,
                           split=".sign") +
                facet_grid(.~.sign) +
                labs(title=NULL)
            save_plot(p_dot, paste0(output_prefix, "_GSEA_Dot_", cid, "_mqc"),
                     12.5, 10.5)

            # Enrichment map
            tryCatch({
                p_emap <- emapplot(gsea_out,
                                  showCategory=GSEA_SHOW_CATEGORIES_EMAP,
                                  cex.params=list(category_label=0.7)) +
                    labs(title=NULL)
                save_plot(p_emap, paste0(output_prefix, "_GSEA_Emap_", cid, "_mqc"),
                         7, 6)
            }, error=function(e) {
                cat(paste0("    WARNING: Emap plot failed: ", e$message, "\n"))
            })

            # Running enrichment score plot (top 3 pathways)
            tryCatch({
                p_running <- gseaplot2(
                    gsea_out,
                    geneSetID = 1:min(GSEA_SHOW_CATEGORIES_RUNNING, nrow(gsea_out)),
                    pvalue_table = TRUE,
                    base_size = 11
                )
                save_plot(p_running, paste0(output_prefix, "_GSEA_Running_", cid, "_mqc"),
                         10, 8)
            }, error=function(e) {
                cat(paste0("    WARNING: Running score plot failed: ", e$message, "\n"))
            })

            # Network plot (cnetplot)
            tryCatch({
                p_cnet <- cnetplot(
                    gsea_out,
                    showCategory = GSEA_SHOW_CATEGORIES_CNET,
                    foldChange = gene_list,
                    colorEdge = TRUE,
                    cex.params = list(category_label = 0.6, gene_label = 0.5)
                ) +
                labs(title = NULL)
                save_plot(p_cnet, paste0(output_prefix, "_GSEA_Network_", cid, "_mqc"),
                         10, 10)
            }, error=function(e) {
                cat(paste0("    WARNING: Network plot failed: ", e$message, "\n"))
            })

        } else {
            cat("    WARNING: No significant GSEA results\n")
        }

        processed_count <- processed_count + 1

    }, error = function(e) {
        cat(paste0("    ERROR: Failed to process contrast ", cid, ": ",
                  e$message, "\n"))
        failed_count <- failed_count + 1
    })
}

cat(paste0("  Summary: ", processed_count, " contrasts processed successfully, ",
          failed_count, " failed\n"))
# ==============================================================================
# 6. GENERATE GEMINI PROMPT FILE
# ==============================================================================
cat("LOG [6/6]: Generating Gemini Prompt Summary...\n")

if (length(prompt_summary) > 0) {
    prompt_text <- paste0(
        "Loaded Results from Pipeline:\n\n",
        paste(unlist(prompt_summary), collapse="\n")
    )
    writeLines(prompt_text, paste0(dirname(output_prefix), "/GEMINI_PROMPTS.txt"))
    cat("  SUCCESS: Gemini prompts file created\n")
} else {
    cat("  WARNING: No GSEA results available for prompt generation\n")
}

cat("\n=== PIPELINE COMPLETE ===\n")
cat(paste0("Total plots generated for ", processed_count, " contrasts\n"))
cat("Plot types per contrast:\n")
cat("  - Volcano plot\n")
cat("  - MA plot\n")
cat("  - Heatmap of top DE genes\n")
cat("  - Boxplots of top genes\n")
cat("  - GSEA dotplot\n")
cat("  - GSEA enrichment map\n")
cat("  - GSEA running score plot\n")
cat("  - GSEA gene-pathway network\n")
cat("\nGlobal QC plots:\n")
cat("  - Phylogenetic tree\n")
cat("  - PCA with subtype annotations\n")
cat("  - Sample correlation heatmap\n")
cat("  - Subtype signature heatmap\n")

# ADDED: Session info for reproducibility
writeLines(capture.output(sessionInfo()), paste0(dirname(output_prefix), "/sessionInfo.txt"))

