#!/usr/bin/env Rscript
# run_pathways_drugs_v4.R
# FIXED: Replaced S4 accessor '@names' with standard 'names()' to prevent crash.

suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(ape); library(ggrepel)
    library(EnsDb.Hsapiens.v86); library(clusterProfiler); library(enrichplot)
    library(GSVA); library(GSEABase); library(ComplexHeatmap); library(circlize)
    library(data.table); library(tidyr); library(stringr); library(igraph); library(ggraph)
})

# Global reproducibility
set.seed(12345)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
PADJ_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 1.0
STRING_SCORE_CUT <- 700
TOP_HUBS_N <- 15

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

# 1. ROBUST MAPPING
map_genes_to_symbols <- function(gene_ids, db = EnsDb.Hsapiens.v86) {
    clean_ids <- sub("\\..*", "", gene_ids)
    # Check if they are already symbols (heuristic)
    if(mean(grepl("^ENSG", clean_ids)) < 0.1) return(clean_ids)
    
    symbols <- mapIds(db, keys=clean_ids, column="SYMBOL", keytype="GENEID", multiVals="first")
    ifelse(is.na(symbols), clean_ids, symbols)
}

# 2. SAVE PLOT
save_mqc <- function(plot_obj, filename_base, w=10, h=8) {
    tryCatch({
        if (inherits(plot_obj, "Heatmap")) {
            pdf(paste0(filename_base, "_mqc.pdf"), width=w, height=h); draw(plot_obj); dev.off()
            png(paste0(filename_base, "_mqc.png"), width=w, height=h, units="in", res=300); draw(plot_obj); dev.off()
        } else {
            ggsave(paste0(filename_base, "_mqc.pdf"), plot_obj, width=w, height=h)
            ggsave(paste0(filename_base, "_mqc.png"), plot_obj, width=w, height=h, dpi=300, bg="white")
        }
    }, error = function(e) { cat("ERROR saving plot:", e$message, "\n") })
}

# 3. REPORTING
report_content <- character()
add_text <- function(header, content) {
    report_content <<- c(report_content, paste0("\n## ", header, "\n\n", content, "\n"))
}

# ==============================================================================
# SETUP
# ==============================================================================
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<5) stop("Usage: script.R <vst> <results_dir> <gmt_dir> <string_dir> <out_prefix>")

vst_file <- args[1]; results_dir <- args[2]; gmt_dir <- args[3]; string_dir <- args[4]; out_prefix <- args[5]

# Init Report
report_content <<- c(paste0("# Integrated Transcriptomic Analysis\n**Date:** ", Sys.Date(), "\n"))

# ==============================================================================
# 1. LOAD STRING DATABASE (Local)
# ==============================================================================
cat("LOG: Loading STRING Database...\n")
link_f <- list.files(string_dir, pattern="protein.links.*.txt.gz", full.names=TRUE)[1]
info_f <- list.files(string_dir, pattern="protein.info.*.txt.gz", full.names=TRUE)[1]

# Load Map
string_map <- fread(info_f, select=c(1, 2))
colnames(string_map) <- c("id", "symbol")
sym2string <- string_map$id; names(sym2string) <- string_map$symbol
string2sym <- string_map$symbol; names(string2sym) <- string_map$id

# Load Network
string_net <- fread(link_f)
if(ncol(string_net) >= 3) colnames(string_net)[1:3] <- c("protein1", "protein2", "combined_score")
string_net <- string_net[combined_score >= STRING_SCORE_CUT]

# ==============================================================================
# 2. DATA LOADING & CLEANING
# ==============================================================================
cat("LOG: Loading VST Data...\n")
mat_vst <- as.matrix(read.table(vst_file, header=TRUE, row.names=1, check.names=FALSE))

# Map symbols using the robust function
mapped_syms <- map_genes_to_symbols(rownames(mat_vst))
mat_vst_sym <- mat_vst
rownames(mat_vst_sym) <- mapped_syms

# Handle duplicates by averaging
mat_vst_sym <- as.data.frame(mat_vst_sym) %>%
    tibble::rownames_to_column("symbol") %>%
    dplyr::filter(!is.na(symbol)) %>%
    group_by(symbol) %>%
    summarise(across(everything(), mean)) %>%
    tibble::column_to_rownames("symbol") %>%
    as.matrix()

cat(paste0("  Mapped to ", nrow(mat_vst_sym), " unique gene symbols\n"))

# ==============================================================================
# 3. CONTRAST LOOP
# ==============================================================================
contrasts <- list.files(file.path(results_dir, "tables/differential"), pattern=".results.tsv", full.names=TRUE)

for(f in contrasts) {
    cid <- sub(".deseq2.results.tsv", "", basename(f))
    cat(paste0("\nLOG: Processing Contrast: ", cid, "...\n"))
    
    # Load Results
    res_df <- read.table(f, header=TRUE, sep="\t", quote="")
    
    # Standardize ID column
    if (grepl("^[0-9]+$", rownames(res_df)[1])) {
        rownames(res_df) <- res_df$gene_id # Fallback if rownames are indices
    }
    
    # Add Symbols if missing (Using robust function)
    if(!"symbol" %in% colnames(res_df)) {
        res_df$symbol <- map_genes_to_symbols(rownames(res_df))
    }
    
    # PREPARE GSEA LIST
    gene_list <- res_df %>%
        dplyr::filter(!is.na(log2FoldChange), !is.na(symbol), is.finite(log2FoldChange)) %>%
        distinct(symbol, .keep_all = TRUE) %>%
        arrange(desc(log2FoldChange)) %>%
        pull(log2FoldChange, name = symbol)
    
    sig_genes <- res_df %>% dplyr::filter(padj < PADJ_CUTOFF, abs(log2FoldChange) > LOG2FC_CUTOFF) %>% pull(symbol)
    sig_genes <- intersect(sig_genes, rownames(mat_vst_sym))
    
    add_text(cid, paste0("Differential analysis identified ", length(sig_genes), " significant genes."))

    # --------------------------------------------------------------------------
    # A. MULTI-DATABASE GSEA LOOP
    # --------------------------------------------------------------------------
    gmts <- list.files(gmt_dir, pattern=".gmt", full.names=TRUE)
    
    for(gmt_path in gmts) {
        db_name <- tools::file_path_sans_ext(basename(gmt_path))
        
        # Skip special files for later
        if(grepl("dsigdb", db_name) || grepl("string", db_name)) next
        
        cat(paste0("  > Running GSEA: ", db_name, "\n"))
        
        # Load GMT
        gmt_data <- tryCatch(read.gmt(gmt_path), error=function(e) NULL)
        if(is.null(gmt_data)) next
        
        # *** FIX IS HERE: Use names(gene_list) instead of gene_list@names ***
        if(length(intersect(names(gene_list), gmt_data$gene)) < 5) {
            cat("    [SKIP] Low overlap with gene list.\n")
            next
        }

        # Run GSEA
        gsea_out <- tryCatch({
            GSEA(gene_list, TERM2GENE = gmt_data, pvalueCutoff = 1, verbose = FALSE, eps = 1e-50)
        }, error = function(e) { return(NULL) })

        if(!is.null(gsea_out) && nrow(gsea_out) > 0) {
            # 1. Dotplot
            p_dot <- dotplot(gsea_out, showCategory=15, split=".sign") + facet_grid(.~.sign) +
                     ggtitle(paste0(db_name, ": ", cid)) + theme(axis.text.y = element_text(size=7))
            save_mqc(p_dot, paste0(out_prefix, "_", cid, "_GSEA_Dot_", db_name))
            
            # 2. Enrichment Map
            tryCatch({
                gsea_sim <- pairwise_termsim(gsea_out)
                p_emap <- emapplot(gsea_sim, showCategory=20, cex.params=list(category_label=0.6))
                save_mqc(p_emap, paste0(out_prefix, "_", cid, "_GSEA_Emap_", db_name))
            }, error=function(e){})
            
            # CSV Export
            write.csv(gsea_out@result, paste0(out_prefix, "_", cid, "_", db_name, ".csv"))
            
            # Report Text
            top <- gsea_out@result %>% arrange(p.adjust) %>% head(1)
            add_text(paste0("Enrichment (", db_name, ")"), 
                     paste0("Top pathway: **", top$ID, "** (NES=", round(top$NES, 2), ")."))
        }
    }

    # --------------------------------------------------------------------------
    # B. DRUG REPURPOSING (DSigDB)
    # --------------------------------------------------------------------------
    cat("  > Drug Discovery (DSigDB)...\n")
    dsig_path <- list.files(gmt_dir, pattern="dsigdb", full.names=TRUE)[1]
    
    if(!is.na(dsig_path)) {
        drug_gmt <- read.gmt(dsig_path)
        
        drug_gsea <- tryCatch({
            GSEA(gene_list, TERM2GENE = drug_gmt, pvalueCutoff = 1, verbose = FALSE, eps = 1e-50)
        }, error=function(e) NULL)
        
        if(!is.null(drug_gsea)) {
            cands <- drug_gsea@result %>% dplyr::filter(NES < -1.0, p.adjust < 0.05) %>% arrange(p.adjust)
            
            if(nrow(cands) > 0) {
                p <- dotplot(drug_gsea, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("DSigDB Drug Candidates")
                save_mqc(p, paste0(out_prefix, "_", cid, "_Drug_Discovery"))
                
                top_d <- cands$ID[1]
                add_text("Drug Candidates", 
                         paste0("Top candidate: **", top_d, "** (NES=", round(cands$NES[1], 2), ")."))
            }
        }
    }

    # --------------------------------------------------------------------------
    # C. STRING PPI NETWORKS
    # --------------------------------------------------------------------------
    cat("  > PPI Network...\n")
    mapped_ids <- sym2string[sig_genes]
    mapped_ids <- mapped_ids[!is.na(mapped_ids)]
    
    if(length(mapped_ids) > 10) {
        sub_net <- string_net[protein1 %in% mapped_ids & protein2 %in% mapped_ids]
        
        if(nrow(sub_net) > 0) {
            g <- graph_from_data_frame(sub_net, directed=FALSE)
            V(g)$string_id <- V(g)$name
            V(g)$name <- string2sym[V(g)$string_id]
            
            deg <- degree(g)
            hubs <- names(sort(deg, decreasing=TRUE)[1:min(TOP_HUBS_N, length(deg))])
            
            comps <- components(g)
            g_main <- induced_subgraph(g, names(comps$membership[comps$membership == which.max(comps$csize)]))
            V(g_main)$type <- ifelse(V(g_main)$name %in% hubs, "Hub", "Node")
            
            p_net <- ggraph(g_main, layout="fr") + 
                     geom_edge_link(alpha=0.2, color="grey70") + 
                     geom_node_point(aes(color=type, size=type)) + 
                     scale_color_manual(values=c("Hub"="#E41A1C", "Node"="#377EB8")) +
                     scale_size_manual(values=c("Hub"=5, "Node"=2)) +
                     geom_node_text(aes(label=ifelse(type=="Hub", name, "")), repel=TRUE, fontface="bold", bg.color="white") +
                     theme_void() + ggtitle(paste0("STRING PPI: ", cid))
            
            save_mqc(p_net, paste0(out_prefix, "_", cid, "_PPI_Network"))
            add_text("Protein Interaction Network", paste0("Identified Hub Genes: ", paste(hubs, collapse=", ")))
        }
    }
}

# Final HTML Report
md_file <- paste0(dirname(out_prefix), "/ANALYSIS_REPORT_FINAL.md")
writeLines(report_content, md_file)
html_content <- paste0("<div class='mqc-custom-content-section'>", 
                       paste(gsub("## ", "<h3>", gsub("\n", "<br>", report_content)), collapse="\n"), 
                       "</div>")
writeLines(html_content, paste0(dirname(out_prefix), "/Analysis_Narrative_mqc.html"))

cat("\nLOG: Pipeline Complete.\n")
