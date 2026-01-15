#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# THE KITCHEN SINK: VISUALIZATION SUITE (Dictionary-Based Mapping)
# ------------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript plot_kitchen_sink.R <deseq_file> <gmt_file> <output_prefix> <counts_dictionary>")
}

deseq_file <- args[1]
gmt_file   <- args[2]
out_prefix <- args[3]
counts_file <- args[4]

# Create output directory
out_dir <- dirname(out_prefix)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(" [1/8] Loading Libraries...")
suppressPackageStartupMessages({
    library(clusterProfiler)
    library(enrichplot)
    library(data.table)
    library(ggplot2)
    # library(org.Hs.eg.db) REMOVED: We use file-based mapping now
})

# Helper to save plots
generated_images <- list()
save_and_store <- function(plot_obj, suffix, title, w, h) {
    base_name <- paste0(out_prefix, suffix)
    ggsave(paste0(base_name, ".pdf"), plot_obj, width=w, height=h)
    png_filename <- paste0(basename(base_name), ".png")
    ggsave(paste0(base_name, ".png"), plot_obj, width=w, height=h, dpi=300, bg="white")
    
    section_id <- tolower(gsub(" ", "_", title))
    generated_images[[section_id]] <<- list(title = title, file = png_filename)
    message(paste("      Saved:", suffix))
}

# ------------------------------------------------------------------------------
# 1. PREPARE DATA & MAP IDs
# ------------------------------------------------------------------------------
message(paste(" [2/8] Reading Data:", deseq_file))
res <- fread(deseq_file)

# Check for Symbol column
possible_names <- c("gene_name", "GeneName", "symbol", "Symbol", "name")
name_col <- intersect(possible_names, colnames(res))[1]

if (is.na(name_col)) {
    message("LOG: No Gene Symbol column found. Checking for Ensembl IDs in 'gene_id'...")
    
    if ("gene_id" %in% colnames(res)) {
        # Check if they look like Ensembl IDs
        example <- res$gene_id[1]
        # UPDATED: Regex to capture ENSG (Human) or ENS (Rat/General)
        if (grepl("^ENS", example)) { 
            message("LOG: Detected Ensembl IDs. Mapping to Symbols using Counts Dictionary...")
            
            # --- FILE BASED MAPPING ---
            dic <- fread(counts_file, select = c("gene_id", "gene_name"))
            
            # Merge
            res <- merge(res, dic, by="gene_id", all.x=TRUE)
            
            # Filter out failures
            n_before <- nrow(res)
            res <- res[!is.na(gene_name) & gene_name != ""]
            message(paste("LOG: Mapped", nrow(res), "genes out of", n_before))
            
            name_col <- "gene_name"
        } else {
             # Assume the gene_id is actually the symbol
             name_col <- "gene_id"
        }
    } else {
        stop("ERROR: Neither 'gene_name' nor 'gene_id' columns found.")
    }
}

# ------------------------------------------------------------------------------
# 2. CREATE RANKED LIST
# ------------------------------------------------------------------------------
stat_col <- "log2FoldChange"
p_col    <- "pvalue"

if (!stat_col %in% colnames(res)) stop(paste("Column missing:", stat_col))

res$rank_metric <- sign(res[[stat_col]]) * -log10(res[[p_col]])
res <- res[!is.na(rank_metric) & !is.infinite(rank_metric)]

# Handle Duplicates: Take the one with highest magnitude metric
res <- res[order(abs(rank_metric), decreasing = TRUE)]
res <- res[!duplicated(get(name_col))]

# Create named vector
gene_list <- setNames(res$rank_metric, res[[name_col]])
gene_list <- sort(gene_list, decreasing = TRUE)

message(paste("LOG: Top gene:", names(gene_list)[1], "=", gene_list[1]))

# ------------------------------------------------------------------------------
# 3. RUN GSEA
# ------------------------------------------------------------------------------
message(" [3/8] Running GSEA...")
gmt <- read.gmt(gmt_file)

# Check overlap
overlap <- length(intersect(names(gene_list), gmt$gene))
message(paste("LOG: Overlap with GMT:", overlap, "genes"))

if (overlap < 5) stop("ERROR: Almost zero overlap between Data and GMT. Mapping likely failed.")

gsea_res <- GSEA(gene_list, TERM2GENE = gmt, pvalueCutoff = 1.0, minGSSize = 5, maxGSSize = 500, verbose = FALSE, eps=1e-50)

if (is.null(gsea_res) || nrow(gsea_res) == 0) stop("GSEA returned no results.")

# Filter Top 50
top_gsea <- gsea_res
if (nrow(top_gsea) > 50) {
    top_gsea@result <- head(top_gsea@result[order(top_gsea@result$p.adjust), ], 50)
}
top_gsea <- pairwise_termsim(top_gsea)

# ------------------------------------------------------------------------------
# 4. GENERATE PLOTS
# ------------------------------------------------------------------------------
message(" [4/8] Generating Plots...")

# 1. Dotplot
p1 <- dotplot(top_gsea, showCategory=15, split=".sign") + facet_grid(.~.sign) + ggtitle("Activated vs Suppressed Pathways")
save_and_store(p1, "_1_Dotplot", "Dotplot Summary", 12, 9)

# 2. Network
p3 <- emapplot(top_gsea, showCategory = 30, cex_label_category=0.6, layout="nicely") + ggtitle("Pathway Functional Modules")
save_and_store(p3, "_3_EnrichMap", "Enrichment Map", 12, 10)

# 3. Cnet
p5 <- cnetplot(top_gsea, categorySize="pvalue", foldChange=gene_list, showCategory=5, circular=TRUE, colorEdge=TRUE) + ggtitle("Top 5 Pathways & Gene Linkages")
save_and_store(p5, "_5_Cnet", "Circular Network", 12, 12)

# 4. Heatplot
p6 <- heatplot(top_gsea, foldChange=gene_list, showCategory=10) + ggtitle("Gene-Pathway Expression Matrix")
save_and_store(p6, "_6_Heatplot", "Heatplot Matrix", 16, 6)

# ------------------------------------------------------------------------------
# 5. MULTIQC YAML (FIXED FOR RELATIVE PATHS)
# ------------------------------------------------------------------------------
message(" [8/8] Generating MultiQC YAML...")
mqc_file <- paste0(dirname(out_prefix), "/pathway_analysis_mqc.yaml")

# Logic: Since MultiQC usually runs one level up (in the results root), 
# we need to prepend the folder name (e.g., "plots/") to the image src.
folder_name <- basename(dirname(out_prefix))

sink(mqc_file)
cat("id: 'pathway_analysis'\nsection_name: 'Pathway Enrichment'\nplot_type: 'html'\ndata: |\n    <div class='row'>\n")
for (id in names(generated_images)) {
    info <- generated_images[[id]]
    # Dynamic Path Fix:
    rel_path <- paste0(folder_name, "/", info$file)
    cat(paste0("        <div class='col-md-6'><h4>", info$title, "</h4><img src='", rel_path, "' style='width:100%'></div>\n"))
}
cat("    </div>\n")
sink()
message(" [Done]")
