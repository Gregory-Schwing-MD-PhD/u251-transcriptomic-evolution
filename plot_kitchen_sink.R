#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# THE KITCHEN SINK: MASTER VISUALIZATION SUITE (v4 - GSVA Fix, _mqc filenames)
# Combines: GSEA (Dot, Tree, Map, Cnet, UpSet), GSVA (Heatmap), Volcano
# ------------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript plot_kitchen_sink.R <deseq_res_file> <vst_counts_file> <gmt_file> <output_prefix> <counts_dictionary>")
}

deseq_file <- args[1]
vst_file   <- args[2]
gmt_file   <- args[3]
out_prefix <- args[4]
counts_file <- args[5]

# Create output directory
out_dir <- dirname(out_prefix)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(" [1/10] Loading Libraries...")
suppressPackageStartupMessages({
    library(clusterProfiler)
    library(enrichplot)
    library(GSVA)
    library(GSEABase)
    library(ggplot2)
    library(ggnewscale)
    library(EnhancedVolcano)
    library(ComplexHeatmap)
    library(circlize)
    library(stringr)
    library(data.table)
    library(GOSemSim)
    
    # Check for UpSet library
    if (!requireNamespace("ggupset", quietly = TRUE)) {
        message("WARNING: 'ggupset' not found. Skipping UpSet plot.")
    } else {
        library(ggupset)
    }
})

# Helper to save plots
generated_images <- list()
save_and_store <- function(plot_obj, suffix, title, w, h) {
    # Add _mqc so MultiQC recognises these as report assets
    base_name <- paste0(out_prefix, suffix, "_mqc")
    
    # Save PDF
    tryCatch({
        if (inherits(plot_obj, "Heatmap")) {
            pdf(paste0(base_name, ".pdf"), width=w, height=h)
            draw(plot_obj)
            dev.off()
        } else {
            ggsave(paste0(base_name, ".pdf"), plot_obj, width=w, height=h)
        }
    }, error = function(e) message(paste("PDF Error:", e$message)))

    # Save PNG for MultiQC
    png_filename <- paste0(basename(base_name), ".png")
    tryCatch({
        if (inherits(plot_obj, "Heatmap")) {
            png(paste0(base_name, ".png"), width=w, height=h, units="in", res=300)
            draw(plot_obj)
            dev.off()
        } else {
            ggsave(paste0(base_name, ".png"), plot_obj, width=w, height=h, dpi=300, bg="white")
        }
        section_id <- tolower(gsub(" ", "_", title))
        generated_images[[section_id]] <<- list(title = title, file = png_filename)
        message(paste("      Saved:", suffix))
    }, error = function(e) message(paste("PNG Error:", e$message)))
}

# ------------------------------------------------------------------------------
# 2. LOAD & PREPARE DATA
# ------------------------------------------------------------------------------
message(paste(" [2/10] Reading Data..."))

# Load Dictionary for Mapping
dic <- fread(counts_file, select = c("gene_id", "gene_name"))

# A. DESeq2 Results
res <- fread(deseq_file)
# ID Mapping Logic
possible_names <- c("gene_name", "GeneName", "symbol", "Symbol", "name")
name_col <- intersect(possible_names, colnames(res))[1]

if (is.na(name_col)) {
    if ("gene_id" %in% colnames(res)) {
        if (grepl("^ENS", res$gene_id[1])) {
            message("LOG: Mapping Ensembl IDs to Symbols (DESeq2)...")
            res <- merge(res, dic, by="gene_id", all.x=TRUE)
            res <- res[!is.na(gene_name) & gene_name != ""]
            name_col <- "gene_name"
        } else {
             name_col <- "gene_id"
        }
    } else {
        stop("ERROR: No gene ID column found in DESeq2 results.")
    }
}

# B. VST Counts (for GSVA) - FIXED MAPPING
message("LOG: Processing VST Matrix (Ensembl -> Symbol)...")
vst_dt <- fread(vst_file)
colnames(vst_dt)[1] <- "gene_id" # Ensure first col is ID

# Check if mapping is needed (look for ENSG prefix)
if (grepl("^ENS", vst_dt$gene_id[1])) {
    # Merge with dictionary
    vst_dt <- merge(vst_dt, dic, by="gene_id", all.x=TRUE)
    
    # Filter valid symbols
    vst_dt <- vst_dt[!is.na(gene_name) & gene_name != ""]
    
    # Aggregate Duplicates (Calculate Mean VST per Symbol)
    numeric_cols <- setdiff(colnames(vst_dt), c("gene_id", "gene_name"))
    vst_aggr <- vst_dt[, lapply(.SD, mean), by=gene_name, .SDcols=numeric_cols]
    
    # Convert to Matrix
    vst_data <- as.matrix(vst_aggr[, -1])
    rownames(vst_data) <- vst_aggr$gene_name
    message(paste("      Mapped", nrow(vst_data), "unique gene symbols."))
} else {
    # Already symbols?
    vst_data <- as.matrix(vst_dt[, -1])
    rownames(vst_data) <- vst_dt$gene_id
}

# ------------------------------------------------------------------------------
# 3. ENHANCED VOLCANO PLOT
# ------------------------------------------------------------------------------
message(" [3/10] Generating Volcano Plot...")

volcano_args <- list(
    toptable = res,
    lab = res[[name_col]],
    x = 'log2FoldChange',
    y = 'pvalue',
    title = 'Differential Expression',
    subtitle = 'Recurrent vs Primary',
    pCutoff = 1e-05,
    FCcutoff = 1.5,
    pointSize = 2.0,
    labSize = 3.0,
    colAlpha = 0.6,
    legendPosition = 'right',
    legendLabSize = 10,
    legendIconSize = 3.0,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    max.overlaps = 30
)
p_vol <- do.call(EnhancedVolcano, volcano_args)
save_and_store(p_vol, "_1_Volcano", "Volcano Plot", 10, 10)

# ------------------------------------------------------------------------------
# 4. PREPARE GSEA RANKING
# ------------------------------------------------------------------------------
stat_col <- "log2FoldChange"
p_col    <- "pvalue"

res$rank_metric <- sign(res[[stat_col]]) * -log10(res[[p_col]])
res <- res[!is.na(rank_metric) & !is.infinite(rank_metric)]
res <- res[order(abs(rank_metric), decreasing = TRUE)]
res <- res[!duplicated(get(name_col))]

gene_list <- setNames(res$rank_metric, res[[name_col]])
gene_list <- sort(gene_list, decreasing = TRUE)

# ------------------------------------------------------------------------------
# 5. RUN GSEA
# ------------------------------------------------------------------------------
message(" [4/10] Running GSEA...")
gmt <- read.gmt(gmt_file)

overlap <- length(intersect(names(gene_list), gmt$gene))
if (overlap < 5) stop("ERROR: < 5 genes overlap between Data and GMT.")

gsea_res <- GSEA(gene_list, TERM2GENE = gmt, pvalueCutoff = 1.0, minGSSize = 5, maxGSSize = 500, verbose = FALSE, eps=1e-50)

if (is.null(gsea_res) || nrow(gsea_res) == 0) stop("GSEA returned no results.")

# Filter Top 60 for Clustering
top_gsea <- gsea_res
if (nrow(top_gsea) > 60) {
    top_gsea@result <- head(top_gsea@result[order(top_gsea@result$p.adjust), ], 60)
}
top_gsea <- pairwise_termsim(top_gsea)

# ------------------------------------------------------------------------------
# 6. GSEA VISUALIZATIONS
# ------------------------------------------------------------------------------
message(" [5/10] Generating GSEA Plots...")

p_dot <- dotplot(top_gsea, showCategory=15, split=".sign", label_format=40) +
         facet_grid(.~.sign) + ggtitle("Activated vs Suppressed Pathways") +
         theme(axis.text.y = element_text(size=8))
save_and_store(p_dot, "_2_Dotplot", "Dotplot Summary", 12, 10)

tryCatch({
    p_tree <- treeplot(top_gsea, 
                       cluster.params = list(method = "ward.D2"), 
                       showCategory = 30, 
                       fontsize=3) +
              ggtitle("Hierarchical Pathway Clustering")
    save_and_store(p_tree, "_3_Treeplot", "Treeplot Clusters", 14, 10)
}, error = function(e) message(paste("Skipping Treeplot:", e$message)))

p_emap <- emapplot(top_gsea, 
                   showCategory = 30, 
                   cex.params = list(category_label=0.6), 
                   layout.params = list(layout="nicely")) +
          ggtitle("Pathway Functional Modules")
save_and_store(p_emap, "_4_EnrichMap", "Enrichment Map", 12, 10)

p_cnet <- cnetplot(top_gsea, 
                   categorySize="pvalue", 
                   showCategory=5,
                   circular=TRUE, 
                   color.params = list(foldChange=gene_list, edge=TRUE),
                   cex.params = list(category_label=0.7, gene_label=0.5)) +
          ggtitle("Top 5 Pathways & Gene Linkages")
save_and_store(p_cnet, "_5_Cnet", "Gene-Concept Network", 14, 14)

if (requireNamespace("ggupset", quietly = TRUE)) {
    p_upset <- upsetplot(top_gsea, n=10) + ggtitle("Gene Overlap Top 10 Pathways")
    save_and_store(p_upset, "_6_UpSet", "UpSet Intersection", 14, 8)
}

# ------------------------------------------------------------------------------
# 7. RUN GSVA
# ------------------------------------------------------------------------------
message(" [6/10] Running GSVA (Sample-Level)...")

# Ensure gene sets are loaded correctly
gene_sets_list <- split(gmt$gene, gmt$term)

# Run GSVA with method-specific parameter
gsva_res <- gsva(vst_data, gene_sets_list, method="gsva", kcdf="Gaussian", min.sz=10, max.sz=500, verbose=FALSE)

# Select variable pathways
pathway_var <- apply(gsva_res, 1, var)
top_pathways <- names(sort(pathway_var, decreasing=TRUE))[1:30]
gsva_subset <- gsva_res[top_pathways, ]

# ------------------------------------------------------------------------------
# 8. GSVA HEATMAP
# ------------------------------------------------------------------------------
message(" [7/10] Generating GSVA Heatmap...")

ht <- Heatmap(gsva_subset,
        name = "GSVA Score",
        column_title = "Sample-Specific Pathway Activity",
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        col = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
)
save_and_store(ht, "_7_GSVA_Heatmap", "GSVA Heatmap", 12, 10)

# ------------------------------------------------------------------------------
# 9. MULTIQC YAML
# ------------------------------------------------------------------------------
message(" [8/10] Generating MultiQC YAML...")
mqc_file <- paste0(dirname(out_prefix), "/pathway_analysis_mqc.yaml")
folder_name <- basename(dirname(out_prefix))

sink(mqc_file)
cat("id: 'pathway_analysis'\nsection_name: 'Pathway Enrichment'\nplot_type: 'html'\ndata: |\n    <div class='row'>\n")
for (id in sort(names(generated_images))) {
    info <- generated_images[[id]]
    rel_path <- paste0(folder_name, "/", info$file)
    cat(paste0("        <div class='col-md-6'><h4>", info$title,
               "</h4><img src='", rel_path, "' style='width:100%'></div>\n"))
}
cat("    </div>\n")
sink()

message(" [Done] All plots generated successfully.")

