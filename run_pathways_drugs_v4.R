#!/usr/bin/env Rscript
# run_pathways_drugs_v4.R
# Final Pipeline: STRING Local PPI, Multi-DB GSEA, DSigDB Drug Discovery
# Output: _mqc.pdf/png plots and _mqc.html report for MultiQC integration.

suppressPackageStartupMessages({
    library(clusterProfiler); library(enrichplot); library(dplyr); library(ggplot2)
    library(org.Hs.eg.db); library(igraph); library(stringr); library(GSEABase)
    library(data.table); library(ggraph); library(DOSE); library(ComplexHeatmap)
})

# Reproducibility
set.seed(12345)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
PADJ_CUT <- 0.05
LFC_CUT <- 1.0
STRING_SCORE_CUT <- 700  # High Confidence Interaction (0-1000)
TOP_HUBS_N <- 15         # Number of Hub Genes to label in network

# ==============================================================================
# REPORTING ENGINE (Markdown + HTML for MultiQC)
# ==============================================================================
# We accumulate text in a vector and write both MD and HTML at the end
report_content <- character()

add_text <- function(header, content) {
    # Store Markdown format
    entry <- paste0("\n## ", header, "\n\n", content, "\n")
    report_content <<- c(report_content, entry)
}

init_report <- function() {
    report_content <<- c(paste0("# Integrated Transcriptomic Analysis\n**Date:** ", Sys.Date(), "\n"))
    add_text("Methods Summary", 
             "Differential expression analysis was performed using DESeq2. Pathway enrichment utilized GSEA against MSigDB Hallmark, C2 (Curated), KEGG, and GO databases. Drug candidates were identified via the DSigDB pharmacogenomic reference. PPI networks were constructed using STRING v12.0 (score > 0.7).")
}

write_reports <- function(prefix) {
    # 1. Write Markdown (Standard)
    md_file <- paste0(dirname(prefix), "/ANALYSIS_REPORT_FINAL.md")
    writeLines(report_content, md_file)
    
    # 2. Write HTML Fragment (For MultiQC Custom Content)
    # We strip the '#' to make it fit nicely in MultiQC
    html_content <- report_content
    html_content <- gsub("## ", "<h3>", html_content)
    html_content <- gsub("\n", "<br>", html_content)
    html_content <- paste0(html_content, "</h3>") # Close tags roughly
    
    # Wrap in a div for MultiQC
    final_html <- c(
        "",
        "<div class='mqc-custom-content-section'>",
        "<h2>Automated Interpretation Report</h2>",
        paste(html_content, collapse="\n"),
        "</div>"
    )
    
    html_file <- paste0(dirname(prefix), "/Analysis_Narrative_mqc.html")
    writeLines(final_html, html_file)
    message("Reports written to: ", md_file, " and ", html_file)
}

# ==============================================================================
# PLOT SAVER (The MQC Suffix Enforcer)
# ==============================================================================
save_mqc <- function(plot_obj, base_name, w=10, h=8) {
    # PDF (Vector)
    ggsave(paste0(base_name, "_mqc.pdf"), plot_obj, width=w, height=h)
    # PNG (Bitmap for quick preview)
    ggsave(paste0(base_name, "_mqc.png"), plot_obj, width=w, height=h, dpi=300, bg="white")
}

# ==============================================================================
# SETUP & DATA LOADING
# ==============================================================================
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<5) stop("Usage: script.R <vst> <results_dir> <gmt_dir> <string_dir> <out_prefix>")

vst_file <- args[1]; results_dir <- args[2]; gmt_dir <- args[3]; string_dir <- args[4]; out_prefix <- args[5]

# Initialize Report
init_report()

# 1. LOAD STRING DATABASE (Fast FileReader)
cat("LOG: Loading STRING Database (Local)...\n")
link_f <- list.files(string_dir, pattern="protein.links.*.txt.gz", full.names=TRUE)[1]
info_f <- list.files(string_dir, pattern="protein.info.*.txt.gz", full.names=TRUE)[1]

if(is.na(link_f) || is.na(info_f)) stop("CRITICAL: STRING database files not found in ", string_dir)

string_map <- fread(info_f, select=c("string_protein_id", "preferred_name"))
colnames(string_map) <- c("id", "symbol")
sym2string <- string_map$id; names(sym2string) <- string_map$symbol
string2sym <- string_map$symbol; names(string2sym) <- string_map$id

cat("LOG: Loading Network Topology...\n")
string_net <- fread(link_f)
string_net <- string_net[combined_score >= STRING_SCORE_CUT]

# 2. LOAD VST MATRIX
cat("LOG: Loading Transcriptome...\n")
vst <- as.matrix(read.table(vst_file, header=TRUE, row.names=1))
try({
    library(EnsDb.Hsapiens.v86)
    if(grepl("^ENSG", rownames(vst)[1])) {
        rownames(vst) <- mapIds(EnsDb.Hsapiens.v86, keys=sub("\\..*","",rownames(vst)), column="SYMBOL", keytype="GENEID", multiVals="first")
    }
}, silent=TRUE)
vst <- as.data.frame(vst) %>% tibble::rownames_to_column("g") %>% 
       filter(!is.na(g)) %>% group_by(g) %>% summarize(across(everything(),mean)) %>% 
       tibble::column_to_rownames("g") %>% as.matrix()

# ==============================================================================
# MAIN ANALYSIS LOOP
# ==============================================================================
contrasts <- list.files(file.path(results_dir, "tables/differential"), pattern=".results.tsv", full.names=TRUE)

for(f in contrasts) {
    cid <- sub(".deseq2.results.tsv", "", basename(f))
    cat(paste0("\nLOG: Processing Contrast: ", cid, "...\n"))
    
    res <- read.table(f, header=TRUE)
    if(!"symbol" %in% colnames(res)) {
         try({res$symbol <- mapIds(EnsDb.Hsapiens.v86, keys=sub("\\..*","",rownames(res)), column="SYMBOL", keytype="GENEID", multiVals="first")}, silent=TRUE)
    }
    
    gene_list <- res %>% filter(!is.na(log2FoldChange), !is.na(symbol)) %>% 
                 arrange(desc(log2FoldChange)) %>% pull(log2FoldChange, name=symbol)
    
    sig_genes <- res %>% filter(padj < PADJ_CUT, abs(log2FoldChange) > LFC_CUT) %>% pull(symbol)
    sig_genes <- sig_genes[!is.na(sig_genes)]
    
    add_text(paste0("Analysis: ", cid), paste0("Differential expression identified ", length(sig_genes), " significant genes."))

    # --------------------------------------------------------------------------
    # A. PATHWAY ENRICHMENT (Multi-DB GSEA)
    # --------------------------------------------------------------------------
    gmts <- list.files(gmt_dir, pattern=".gmt", full.names=TRUE)
    
    for(gmt in gmts) {
        db_clean <- tools::file_path_sans_ext(basename(gmt))
        if(grepl("dsigdb", db_clean) || grepl("string", db_clean)) next
        
        cat(paste0("  > GSEA: ", db_clean, "\n"))
        term2gene <- read.gmt(gmt)
        ems <- GSEA(gene_list, TERM2GENE=term2gene, pvalueCutoff=1, verbose=FALSE, eps=1e-50)
        
        if(nrow(ems) > 0) {
            p <- dotplot(ems, showCategory=12, split=".sign") + facet_grid(.~.sign) + 
                 ggtitle(paste0(db_clean, ": ", cid)) + theme(axis.text.y = element_text(size=6))
            
            # SAVE WITH _MQC SUFFIX
            save_mqc(p, paste0(out_prefix, "_", cid, "_", db_clean))
            write.csv(ems@result, paste0(out_prefix, "_", cid, "_", db_clean, ".csv"))
            
            top <- ems@result %>% arrange(p.adjust) %>% head(1)
            add_text(paste0("Enrichment (", db_clean, ")"),
                     paste0("Top pathway: **", top$ID, "** (NES=", round(top$NES, 2), ")."))
        }
    }

    # --------------------------------------------------------------------------
    # B. DRUG REPURPOSING (DSigDB)
    # --------------------------------------------------------------------------
    cat("  > Drug Discovery (DSigDB)...\n")
    dsig_file <- list.files(gmt_dir, pattern="dsigdb", full.names=TRUE)[1]
    
    if(!is.na(dsig_file)) {
        drug_ems <- GSEA(gene_list, TERM2GENE=read.gmt(dsig_file), pvalueCutoff=1, verbose=FALSE, eps=1e-50)
        cands <- drug_ems@result %>% filter(NES < -1.5, p.adjust < 0.05) %>% arrange(p.adjust)
        
        if(nrow(cands) > 0) {
            top_d <- cands$ID[1]
            p <- dotplot(drug_ems, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("Drug Candidates (DSigDB)")
            
            # SAVE WITH _MQC SUFFIX
            save_mqc(p, paste0(out_prefix, "_", cid, "_Drug_GSEA"))
            
            add_text("Pharmacogenomic Candidates",
                     paste0("Analysis against the DSigDB identified **", top_d, "** as a primary candidate. ",
                            "The tumor expression profile correlates negatively with this drug's signature (NES=", round(cands$NES[1], 2), ")."))
        }
    }

    # --------------------------------------------------------------------------
    # C. STRING PROTEIN-PROTEIN INTERACTION NETWORK
    # --------------------------------------------------------------------------
    cat("  > Constructing STRING PPI Network...\n")
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
                     scale_size_manual(values=c("Hub"=6, "Node"=2.5)) +
                     geom_node_text(aes(label=ifelse(type=="Hub", name, "")), 
                                    repel=TRUE, fontface="bold", color="black", bg.color="white", bg.r=0.15) +
                     theme_void() + 
                     ggtitle(paste0("STRING PPI Network: ", cid))
            
            # SAVE WITH _MQC SUFFIX
            save_mqc(p_net, paste0(out_prefix, "_", cid, "_STRING_PPI"))
            
            add_text("PPI Network Topology", 
                        paste0("A high-confidence PPI network was constructed. ",
                               "Topological analysis identified **", paste(hubs, collapse=", "), "** as the central regulatory hubs."))
        }
    }
}

# 3. WRITE FINAL REPORTS
write_reports(out_prefix)
cat("\nLOG: Pipeline Complete. HTML Report and MQC Plots generated.\n")
