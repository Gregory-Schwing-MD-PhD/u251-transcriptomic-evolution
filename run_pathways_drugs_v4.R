#!/usr/bin/env Rscript
# run_pathways_drugs_v4.R (v17 - Config Match & Size Filters)
# UPDATES:
# - Added minGSSize=15 and maxGSSize=500 to match your previous config.
# - Ranks by 'stat' (Wald) -> This is the superior version of "Diff_of_Classes" for RNA-seq.

suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(ape); library(ggrepel)
    library(EnsDb.Hsapiens.v86); library(clusterProfiler); library(enrichplot)
    library(GSVA); library(GSEABase); library(ComplexHeatmap); library(circlize)
    library(data.table); library(tidyr); library(stringr); library(igraph); library(ggraph)
})

# Global Set Seed (Plus seed=TRUE in GSEA calls)
set.seed(12345)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
PADJ_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 1.0
STRING_SCORE_CUT <- 400
TOP_HUBS_N <- 15
PLOT_W <- 16
PLOT_H <- 16
SHOW_CAT_N <- 30
TEXT_REPORT_N <- 10

# GSEA SPECIFIC SETTINGS (Matched to your request)
GSEA_MIN_SIZE <- 15   # Filter out tiny sets (<15 genes)
GSEA_MAX_SIZE <- 500  # Filter out huge sets (>500 genes)

# ==============================================================================
# REPORTING ENGINE
# ==============================================================================
html_buffer <- character()

init_html <- function() {
    style <- "
    <style>
        .mqc-custom-content-section { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; color: #333; }
        .contrast-block { background: #fff; border: 1px solid #ddd; border-radius: 5px; margin-bottom: 40px; padding: 25px; box-shadow: 0 4px 10px rgba(0,0,0,0.05); }
        .contrast-title { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; margin-top: 0; font-size: 24px; }
        .stat-line { font-size: 15px; color: #7f8c8d; margin-bottom: 20px; font-style: italic; }
        .db-header { font-weight: bold; font-size: 16px; color: #fff; background-color: #34495e; padding: 8px 15px; border-radius: 4px 4px 0 0; margin-top: 20px; }
        .pathway-table { width: 100%; border-collapse: collapse; font-size: 12px; border: 1px solid #eee; }
        .pathway-table td, .pathway-table th { padding: 5px 8px; border-bottom: 1px solid #eee; text-align: left; }
        .pathway-table th { background-color: #f8f9fa; color: #555; }
        .nes-pos { color: #c0392b; font-weight: bold; }
        .nes-neg { color: #2980b9; font-weight: bold; }
        .hub-box { padding:10px; background:#f9f9f9; border:1px solid #eee; font-family: monospace; color: #d35400; font-size: 13px;}
        .legend-box { background: #eef6fc; padding: 15px; border-radius: 5px; border: 1px solid #dce4ec; margin-bottom: 30px; }
        .sig-green { color: #27ae60; font-weight: bold; }
        .sig-orange { color: #e67e22; font-style: italic; }
    </style>
    <div class='mqc-custom-content-section'>
    <h1>Integrated Transcriptomic Analysis Report</h1>
    <p><strong>Date:</strong> "
    
    legend <- "
    <div class='legend-box'>
        <div><strong>Interpretation Guide:</strong></div>
        <ul style='margin-top:5px; padding-left:20px;'>
            <li><strong>NES (Score):</strong> <span class='nes-pos'>Red</span> = Activated, <span class='nes-neg'>Blue</span> = Suppressed.</li>
            <li><strong>FDR (Significance):</strong> <span class='sig-green'>Green</span> = Significant (p &lt; 0.05), <span class='sig-orange'>Orange</span> = Exploratory (p &gt; 0.05).</li>
        </ul>
    </div>"
    
    html_buffer <<- c(html_buffer, paste0(style, Sys.Date(), "</p>", legend))
}

add_contrast_header <- function(cid, n_genes, metric) {
    block <- paste0(
        "<div class='contrast-block'>",
        "<h2 class='contrast-title'>", cid, "</h2>",
        "<div class='stat-line'>Differential expression identified <strong>", n_genes, "</strong> significant genes (FDR < ", PADJ_CUTOFF, ").</div>",
        "<div class='stat-line' style='font-size:12px; color:#555;'>GSEA Ranking Metric: <strong>", metric, "</strong> (Selected for small sample size)</div>"
    )
    html_buffer <<- c(html_buffer, block)
}

add_table_from_df <- function(db, df, is_drug=FALSE) {
    header_style <- if(is_drug) "background-color: #27ae60;" else ""
    
    if(is.null(df) || nrow(df) == 0) {
        html_buffer <<- c(html_buffer, paste0("<div class='db-header' style='", header_style, "'>", db, "</div><div style='padding:10px; color:#999;'>No results found.</div>"))
        return()
    }
    
    rows <- apply(df, 1, function(r) {
        nes <- as.numeric(r['NES'])
        fdr <- as.numeric(r['p.adjust'])
        nes_class <- ifelse(nes > 0, "nes-pos", "nes-neg")
        fdr_class <- ifelse(fdr < 0.05, "sig-green", "sig-orange")
        paste0("<tr><td>", r['ID'], "</td>",
               "<td><span class='", nes_class, "'>", round(nes, 2), "</span></td>",
               "<td class='", fdr_class, "'>", formatC(fdr, format="e", digits=2), "</td></tr>")
    })
    
    table_html <- paste0("<div class='db-header' style='", header_style, "'>", db, "</div>",
                         "<table class='pathway-table'><tr><th>ID</th><th>NES</th><th>FDR</th></tr>",
                         paste(rows, collapse=""), "</table>")
    html_buffer <<- c(html_buffer, table_html)
}

add_hubs <- function(hubs) {
    header <- "<div class='db-header' style='background-color: #8e44ad;'>PPI Hub Genes</div>"
    content <- if(is.null(hubs)) "No hubs." else paste0("<div class='hub-box'>", paste(hubs, collapse=", "), "</div>")
    html_buffer <<- c(html_buffer, paste0(header, content))
}

close_block <- function() { html_buffer <<- c(html_buffer, "</div>") }
finish_html <- function(prefix) { html_buffer <<- c(html_buffer, "</div>"); writeLines(html_buffer, paste0(dirname(prefix), "/Analysis_Narrative_mqc.html")) }

# ==============================================================================
# UTILITY
# ==============================================================================
map_genes_to_symbols <- function(gene_ids, db = EnsDb.Hsapiens.v86) {
    clean_ids <- sub("\\..*", "", gene_ids)
    if(mean(grepl("^ENSG", clean_ids)) < 0.1) return(clean_ids)
    symbols <- mapIds(db, keys=clean_ids, column="SYMBOL", keytype="GENEID", multiVals="first")
    ifelse(is.na(symbols), clean_ids, symbols)
}

save_mqc <- function(plot_obj, filename_base, w=PLOT_W, h=PLOT_H) {
    tryCatch({
        ggsave(paste0(filename_base, "_mqc.pdf"), plot_obj, width=w, height=h)
        ggsave(paste0(filename_base, "_mqc.png"), plot_obj, width=w, height=h, dpi=300, bg="white")
    }, error = function(e) { cat("ERROR saving plot:", e$message, "\n") })
}

# ==============================================================================
# MAIN
# ==============================================================================
args <- commandArgs(trailingOnly=TRUE)
vst_file <- args[1]; results_dir <- args[2]; gmt_dir <- args[3]; string_dir <- args[4]; out_prefix <- args[5]

init_html()

cat("LOG: Loading STRING...\n")
link_f <- list.files(string_dir, pattern="protein.links.*.txt.gz", full.names=TRUE)[1]
info_f <- list.files(string_dir, pattern="protein.info.*.txt.gz", full.names=TRUE)[1]
string_map <- fread(info_f, select=c(1, 2)); colnames(string_map) <- c("id", "symbol")
sym2string <- string_map$id; names(sym2string) <- string_map$symbol
string2sym <- string_map$symbol; names(string2sym) <- string_map$id
string_net <- fread(link_f)
if(ncol(string_net) >= 3) colnames(string_net)[1:3] <- c("protein1", "protein2", "combined_score")
string_net <- string_net[combined_score >= STRING_SCORE_CUT]

cat("LOG: Loading VST...\n")
mat_vst <- as.matrix(read.table(vst_file, header=TRUE, row.names=1, check.names=FALSE))
rownames(mat_vst) <- map_genes_to_symbols(rownames(mat_vst))
mat_vst_sym <- as.data.frame(mat_vst) %>% tibble::rownames_to_column("symbol") %>%
    dplyr::filter(!is.na(symbol)) %>% group_by(symbol) %>% summarise(across(everything(), mean)) %>%
    tibble::column_to_rownames("symbol") %>% as.matrix()

contrasts <- list.files(file.path(results_dir, "tables/differential"), pattern=".results.tsv", full.names=TRUE)

for(f in contrasts) {
    cid <- sub(".deseq2.results.tsv", "", basename(f))
    cat(paste0("\nLOG: Contrast: ", cid, "...\n"))
    
    res_df <- read.table(f, header=TRUE, sep="\t", quote="")
    if (grepl("^[0-9]+$", rownames(res_df)[1])) rownames(res_df) <- res_df$gene_id 
    if(!"symbol" %in% colnames(res_df)) res_df$symbol <- map_genes_to_symbols(rownames(res_df))
    
    # --------------------------------------------------------------------------
    # SMART RANKING (Optimized for Small N)
    # --------------------------------------------------------------------------
    metric_name <- "log2FoldChange"
    if("stat" %in% colnames(res_df)) {
        metric_name <- "Wald Statistic (stat)"
        res_df$rank_metric <- res_df$stat
    } else if ("pvalue" %in% colnames(res_df)) {
        metric_name <- "Signed P-value"
        res_df$rank_metric <- sign(res_df$log2FoldChange) * -log10(res_df$pvalue)
    } else {
        res_df$rank_metric <- res_df$log2FoldChange
    }
    cat(paste0("  > Using Ranking Metric: ", metric_name, "\n"))

    # Create Ranked List
    gene_list <- res_df %>%
        dplyr::filter(!is.na(rank_metric), !is.na(symbol), is.finite(rank_metric)) %>%
        distinct(symbol, .keep_all = TRUE) %>%
        arrange(desc(rank_metric)) %>%
        pull(rank_metric, name = symbol)
    
    sig_genes <- res_df %>% dplyr::filter(padj < PADJ_CUTOFF, abs(log2FoldChange) > LOG2FC_CUTOFF) %>% pull(symbol)
    sig_genes <- intersect(sig_genes, rownames(mat_vst_sym))
    
    add_contrast_header(cid, length(sig_genes), metric_name)

    # A. PATHWAYS
    gmts <- list.files(gmt_dir, pattern=".gmt", full.names=TRUE)
    for(gmt_path in gmts) {
        db_name <- tools::file_path_sans_ext(basename(gmt_path))
        if(grepl("dsigdb", db_name) || grepl("string", db_name)) next
        cat(paste0("  > Running GSEA: ", db_name, "\n"))
        
        gmt_data <- tryCatch(read.gmt(gmt_path), error=function(e) NULL)
        if(is.null(gmt_data) || length(intersect(names(gene_list), gmt_data$gene)) < 5) next
        
        # *** FIX: Added minGSSize/maxGSSize and seed ***
        gsea_out <- tryCatch(GSEA(gene_list, TERM2GENE=gmt_data, pvalueCutoff=1, 
                                  minGSSize = GSEA_MIN_SIZE, maxGSSize = GSEA_MAX_SIZE,
                                  verbose=FALSE, eps=1e-50, seed=TRUE), error=function(e) NULL)
        
        if(!is.null(gsea_out) && nrow(gsea_out) > 0) {
            p_dot <- dotplot(gsea_out, showCategory=SHOW_CAT_N, split=".sign") + facet_grid(.~.sign) +
                     ggtitle(paste0(db_name, ": ", cid)) + theme(axis.text.y = element_text(size=7))
            save_mqc(p_dot, paste0(out_prefix, "_", cid, "_GSEA_Dot_", db_name))
            write.csv(gsea_out@result, paste0(out_prefix, "_", cid, "_", db_name, ".csv"))
            
            res <- gsea_out@result
            top_pos <- res %>% filter(NES > 0) %>% arrange(desc(NES)) %>% head(SHOW_CAT_N)
            top_neg <- res %>% filter(NES < 0) %>% arrange(NES) %>% head(SHOW_CAT_N)
            add_table_from_df(db_name, bind_rows(top_pos, top_neg), is_drug=FALSE)
        } else {
            add_table_from_df(db_name, NULL)
        }
    }

    # B. DRUGS
    dsig_path <- list.files(gmt_dir, pattern="dsigdb", full.names=TRUE)[1]
    if(!is.na(dsig_path)) {
        cat(paste0("  > Running Drug Discovery...\n"))
        drug_gmt <- read.gmt(dsig_path)
        # *** FIX: Added Size Filters + Seed ***
        drug_gsea <- tryCatch(GSEA(gene_list, TERM2GENE=drug_gmt, pvalueCutoff=1, 
                                   minGSSize = GSEA_MIN_SIZE, maxGSSize = GSEA_MAX_SIZE,
                                   verbose=FALSE, eps=1e-50, seed=TRUE), error=function(e) NULL)
        
        if(!is.null(drug_gsea)) {
            res <- drug_gsea@result
            # Unified Logic: Top 30 Negative NES
            top_cands <- res %>% filter(NES < 0) %>% arrange(NES) %>% head(SHOW_CAT_N)
            
            if(nrow(top_cands) > 0) {
                p <- dotplot(drug_gsea, showCategory=SHOW_CAT_N, split=".sign") + facet_grid(.~.sign) + ggtitle(paste0("Drug Candidates: ", cid))
                save_mqc(p, paste0(out_prefix, "_", cid, "_Drug_Discovery"))
                add_table_from_df("Therapeutic Candidates (DSigDB)", top_cands, is_drug=TRUE)
            } else {
                add_table_from_df("Therapeutic Candidates (DSigDB)", NULL, is_drug=TRUE)
            }
        }
    }

    # C. PPI
    mapped_ids <- sym2string[sig_genes]
    mapped_ids <- mapped_ids[!is.na(mapped_ids)]
    hub_list <- NULL
    if(length(mapped_ids) > 5) {
        cat(paste0("  > Running PPI Network...\n"))
        sub_net <- string_net[protein1 %in% mapped_ids & protein2 %in% mapped_ids]
        if(nrow(sub_net) > 0) {
            g <- graph_from_data_frame(sub_net, directed=FALSE)
            V(g)$string_id <- V(g)$name
            V(g)$name <- string2sym[V(g)$string_id]
            deg <- degree(g)
            hub_list <- names(sort(deg, decreasing=TRUE)[1:min(TOP_HUBS_N, length(deg))])
            comps <- components(g)
            g_main <- induced_subgraph(g, names(comps$membership[comps$membership == which.max(comps$csize)]))
            V(g_main)$type <- ifelse(V(g_main)$name %in% hub_list, "Hub", "Node")
            p_net <- ggraph(g_main, layout="fr") + 
                     geom_edge_link(alpha=0.2, color="grey70") + 
                     geom_node_point(aes(color=type, size=type)) + 
                     scale_color_manual(values=c("Hub"="#E41A1C", "Node"="#377EB8")) + 
                     scale_size_manual(values=c("Hub"=5, "Node"=2)) +
                     geom_node_text(aes(label=ifelse(type=="Hub", name, "")), repel=TRUE, fontface="bold", bg.color="white") +
                     theme_void() + ggtitle(paste0("PPI: ", cid))
            save_mqc(p_net, paste0(out_prefix, "_", cid, "_PPI_Network"))
        }
    }
    add_hubs(hub_list)
    close_block()
}

finish_html(out_prefix)
cat("\nLOG: Pipeline Complete.\n")
