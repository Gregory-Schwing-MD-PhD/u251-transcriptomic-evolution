#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# EVOLUTIONARY KITCHEN SINK v15.2 (Fixed ID Loading in Step 4)
# ------------------------------------------------------------------------------
# 1. PCA/Tree: Uses 'all.vst.tsv' calculated by the pipeline.
# 2. Volcano:  Uses 'deseq2.results.tsv' calculated by the pipeline.
# 3. GSEA:     Uses Log2FC from the pipeline to generate network plots.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(ape); library(ggrepel)
    library(EnsDb.Hsapiens.v86); library(clusterProfiler); library(enrichplot)
    library(GSVA); library(GSEABase); library(ComplexHeatmap); library(circlize)
    library(EnhancedVolcano); library(data.table); library(tidyr)
    if (requireNamespace("ggupset", quietly = TRUE)) library(ggupset)
})

# --- ARGS ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Usage: script.R <results_dir> <out_prefix> <meta> <contrasts> <gmt>")

results_dir   <- args[1]
output_prefix <- args[2]
meta_file     <- args[3]
contrast_file <- args[4]
gmt_file      <- args[5]

# Output setup
dir.create(dirname(output_prefix), showWarnings = FALSE, recursive = TRUE)
gemini_dir <- paste0(dirname(output_prefix), "/GEMINI_DATA")
dir.create(gemini_dir, showWarnings = FALSE)

# --- SAVE FUNCTION ---
save_plot <- function(plot_obj, filename_base, w=12, h=10) {
    tryCatch({
        if (inherits(plot_obj, "Heatmap")) {
            pdf(paste0(filename_base, ".pdf"), width=w, height=h); draw(plot_obj); dev.off()
            png(paste0(filename_base, ".png"), width=w, height=h, units="in", res=300); draw(plot_obj); dev.off()
        } else {
            ggsave(paste0(filename_base, ".pdf"), plot_obj, width=w, height=h)
            ggsave(paste0(filename_base, ".png"), plot_obj, width=w, height=h, dpi=300, bg="white")
        }
        cat(paste0("SUCCESS: Saved ", basename(filename_base), "\n"))
    }, error = function(e) cat(paste0("ERROR: Failed to save ", basename(filename_base), ": ", e$message, "\n")))
}

# ==============================================================================
# 1. LOAD PRE-COMPUTED VST
# ==============================================================================
cat("LOG [1/5]: Loading Pipeline VST Matrix...\n")
vst_file <- file.path(results_dir, "tables/processed_abundance/all.vst.tsv")

if (!file.exists(vst_file)) stop(paste("FATAL: Pipeline VST file not found at:", vst_file))

# Load VST
mat_vst <- as.matrix(read.table(vst_file, header=TRUE, row.names=1, check.names=FALSE))

# Load Metadata & Sync
meta <- read.csv(meta_file, row.names=1)
common <- intersect(colnames(mat_vst), rownames(meta))
mat_vst <- mat_vst[, common]; meta <- meta[common, , drop=FALSE]

# ==============================================================================
# 2. EVOLUTIONARY TREE
# ==============================================================================
cat("LOG [2/5]: Generating Tree...\n")
top_var <- head(order(apply(mat_vst, 1, var), decreasing=TRUE), 500)
mat_sig <- mat_vst[top_var, ]

phylo_tree <- as.phylo(hclust(dist(t(mat_sig)), method="ward.D2"))
tryCatch({
    if("C2B" %in% phylo_tree$tip.label) rooted_tree <- root(phylo_tree, outgroup="C2B", resolve.root=TRUE) else rooted_tree <- phylo_tree
    pdf(paste0(output_prefix, "_Phylogenetic_Tree_mqc.pdf"), width=10, height=8)
    plot(rooted_tree, main="Evolutionary Reconstruction", type="phylogram", edge.width=2, cex=1.2)
    dev.off()
}, error = function(e) cat("WARN: Tree failed.\n"))

# ==============================================================================
# 3. PCA & SUBTYPING
# ==============================================================================
cat("LOG [3/5]: Generating PCA...\n")
cat("   -> Mapping Ensembl IDs to Gene Symbols for GSVA...\n")

clean_ids <- sub("\\..*", "", rownames(mat_vst))
mapped_syms <- mapIds(EnsDb.Hsapiens.v86, keys=clean_ids, column="SYMBOL", keytype="GENEID", multiVals="first")

mat_vst_sym <- mat_vst
rownames(mat_vst_sym) <- ifelse(is.na(mapped_syms), rownames(mat_vst), mapped_syms)
mat_vst_sym <- mat_vst_sym[!is.na(rownames(mat_vst_sym)), ]
mat_vst_sym <- mat_vst_sym[!duplicated(rownames(mat_vst_sym)), ]

gbm_sigs <- list(
    Mesenchymal = c("CHI3L1", "CD44", "VIM", "RELB", "STAT3", "MET", "TRADD", "MMP9", "TIMP1"),
    Classical   = c("EGFR", "AKT2", "NOTCH3", "JAG1", "CCND2", "F3", "PDGFA", "NES"),
    Proneural   = c("PDGFRA", "IDH1", "OLIG2", "SOX2", "NKX2-2", "OLIG1", "TP53"),
    Neural      = c("NEFL", "GABRA1", "SYT1", "SLC12A5", "MBP", "GABRG2")
)

gsva_sub <- gsva(mat_vst_sym, gbm_sigs, method="gsva", kcdf="Gaussian", verbose=FALSE)

pca <- prcomp(t(mat_sig))
pcaData <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Sample=rownames(meta), Class=meta$Classification)
pcaData$Dominant_Subtype <- apply(gsva_sub, 2, function(x) rownames(gsva_sub)[which.max(x)])
write.csv(pcaData, paste0(gemini_dir, "/pca_coordinates.csv"))

subtype_colors <- c("Mesenchymal"="#E41A1C", "Proneural"="#377EB8", "Classical"="#4DAF4A", "Neural"="#984EA3", "Other/Unclassified"="grey70")
p_pca <- ggplot(pcaData, aes(x=PC1, y=PC2)) +
    geom_point(aes(fill=Dominant_Subtype, shape=Class), size=8, color="black") +
    scale_fill_manual(values=subtype_colors) +
    scale_shape_manual(values=c(21, 24, 22)) +
    theme_bw(base_size=16) +
    labs(title="Evolutionary Trajectory (Pipeline VST)")
save_plot(p_pca, paste0(output_prefix, "_PCA_Publication_mqc"), 10, 8)

# ==============================================================================
# 4. CONTRAST VISUALIZATION
# ==============================================================================
cat("LOG [4/5]: Visualizing Pipeline DE Tables...\n")
contrasts <- read.csv(contrast_file)
gmt <- read.gmt(gmt_file)
prompt_summary <- list()

for(i in 1:nrow(contrasts)) {
    cid <- contrasts$id[i]; target <- contrasts$target[i]; ref <- contrasts$reference[i]
    res_file <- file.path(results_dir, "tables/differential", paste0(cid, ".deseq2.results.tsv"))
    
    cat(sprintf("    -> Processing: %s\n", cid))
    if(!file.exists(res_file)) { cat(sprintf("        WARN: Missing %s\n", res_file)); next }

    # LOAD & CHECK IDS
    # Detect if IDs are in a column or rownames
    res_df <- read.table(res_file, header=TRUE, sep="\t", quote="")
    
    # Check if rownames are just numbers (1, 2, 3...)
    if (grepl("^[0-9]+$", rownames(res_df)[1])) {
        # Try to find the ID column
        id_col_idx <- which(colnames(res_df) %in% c("gene_id", "gene_name", "id", "GeneID"))
        if(length(id_col_idx) > 0) {
            rownames(res_df) <- res_df[, id_col_idx[1]]
        } else {
            # Assume first column is ID
            rownames(res_df) <- res_df[,1]
        }
    }

    # Add symbols if missing
    if(!"symbol" %in% colnames(res_df)) {
        raw_ids <- rownames(res_df)
        clean_ids <- sub("\\..*", "", raw_ids) # Remove version .1
        clean_ids <- gsub('"', '', clean_ids)  # Remove quotes
        
        # DEBUG LOG
        if (i == 1) cat(paste0("       DEBUG: First 5 IDs seen: ", paste(head(clean_ids, 5), collapse=", "), "\n"))
        
        syms <- mapIds(EnsDb.Hsapiens.v86, keys=clean_ids, column="SYMBOL", keytype="GENEID", multiVals="first")
        res_df$symbol <- ifelse(is.na(syms), clean_ids, syms)
    }

    # EXPORT RAW
    write.csv(head(res_df %>% arrange(padj), 50), paste0(gemini_dir, "/DE_Genes_", cid, ".csv"))

    # 1. VOLCANO
    p_vol <- EnhancedVolcano(res_df, lab=res_df$symbol, x='log2FoldChange', y='padj',
        title=paste0("Volcano: ", cid), subtitle=paste0(target, " vs ", ref),
        pCutoff=0.05, FCcutoff=1.0, pointSize=2.0, labSize=3.0,
        drawConnectors=TRUE, max.overlaps=30)
    save_plot(p_vol, paste0(output_prefix, "_Volcano_", cid, "_mqc"), 10, 10)

    # 2. GSEA
    res_df$rank <- res_df$log2FoldChange
    res_df <- res_df[!is.na(res_df$rank) & !is.infinite(res_df$rank),]
    gene_list <- sort(setNames(res_df$rank, res_df$symbol), decreasing=TRUE)

    gsea_out <- GSEA(gene_list, TERM2GENE=gmt, pvalueCutoff=1, verbose=FALSE, eps=1e-50)

    if(!is.null(gsea_out) && nrow(gsea_out) > 0) {
        if(nrow(gsea_out) > 50) gsea_out@result <- head(gsea_out@result[order(gsea_out@result$p.adjust), ], 50)
        write.csv(gsea_out@result[,c("ID", "NES", "p.adjust", "core_enrichment")], paste0(gemini_dir, "/GSEA_Results_", cid, ".csv"))

        top_up <- head(gsea_out@result[gsea_out@result$NES > 0, "ID"], 5)
        top_dn <- head(gsea_out@result[gsea_out@result$NES < 0, "ID"], 5)
        prompt_summary[[cid]] <- paste0("Contrast ", cid, ":\n- Top Up: ", paste(top_up, collapse=", "), "\n- Top Dn: ", paste(top_dn, collapse=", "), "\n")

        gsea_out <- pairwise_termsim(gsea_out)
        save_plot(dotplot(gsea_out, showCategory=15, split=".sign") + facet_grid(.~.sign), paste0(output_prefix, "_GSEA_Dot_", cid, "_mqc"), 12, 16)
        tryCatch({ save_plot(cnetplot(gsea_out, categorySize="pvalue", foldChange=gene_list, circular=TRUE, showCategory=8), paste0(output_prefix, "_GSEA_Cnet_", cid, "_mqc"), 18, 18) }, error=function(e) NULL)
        tryCatch({ save_plot(emapplot(gsea_out, showCategory=20, cex.params=list(category_label=0.7)), paste0(output_prefix, "_GSEA_Emap_", cid, "_mqc"), 16, 14) }, error=function(e) NULL)
    }
}

cat("LOG [5/5]: Generating Prompt...\n")
writeLines(paste0("Loaded Results from Pipeline:\n\n", paste(unlist(prompt_summary), collapse="\n")), paste0(dirname(output_prefix), "/GEMINI_PROMPTS.txt"))
