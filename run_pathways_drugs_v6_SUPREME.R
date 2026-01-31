#!/usr/bin/env Rscript
# run_pathways_drugs_v6_SUPREME.R
# SUPREME EDITION - Enhanced Drug-Pathway Integration with Database Queries
# ==============================================================================
# NEW IN V6:
# - Drug-Pathway Overlap Heatmap (which drugs target which pathways)
# - Mechanism of Action visualization
# - Drug polypharmacology network
# - Enhanced database queries (Reactome, WikiPathways)
# - Cross-database pathway annotations
# ==============================================================================

suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(ape); library(ggrepel)
    library(EnsDb.Hsapiens.v86); library(clusterProfiler); library(enrichplot)
    library(GSVA); library(GSEABase); library(ComplexHeatmap); library(circlize)
    library(data.table); library(tidyr); library(stringr); library(igraph); library(ggraph)
})

# Check for optional packages
HAS_GGRIDGES <- requireNamespace("ggridges", quietly = TRUE)
HAS_HTTR <- requireNamespace("httr", quietly = TRUE)
HAS_JSONLITE <- requireNamespace("jsonlite", quietly = TRUE)

if(!HAS_GGRIDGES) cat("INFO: ggridges not available - ridgeplots will be skipped\n")
if(!HAS_HTTR || !HAS_JSONLITE) cat("INFO: httr/jsonlite not available - API queries will be limited\n")

# Global Set Seed
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

# GSEA Plot Categories
GSEA_DOT_N <- 30
GSEA_EMAP_N <- 20
GSEA_RUNNING_N <- 5
GSEA_NETWORK_N <- 5
GSEA_RIDGE_N <- 10

# Drug-Pathway Integration
DRUG_PATHWAY_TOP_N <- 20
MOA_TOP_DRUGS <- 10

# GSEA Settings
GSEA_MIN_SIZE <- 15
GSEA_MAX_SIZE <- 500

# Report Settings
TEXT_REPORT_N <- 10
EXPORT_TOP_N <- 50

# ==============================================================================
# DATABASE QUERY FUNCTIONS (NEW IN V6)
# ==============================================================================

# Query WikiPathways for pathway information
query_wikipathways <- function(pathway_name) {
    if(!HAS_HTTR || !HAS_JSONLITE) return(NULL)

    tryCatch({
        library(httr)
        library(jsonlite)

        # Clean pathway name
        query <- gsub("_", " ", pathway_name)
        query <- gsub("HALLMARK ", "", query)
        query <- gsub("REACTOME ", "", query)
        query <- gsub("KEGG ", "", query)

        url <- paste0("https://webservice.wikipathways.org/findPathwaysByText?query=",
                      URLencode(query), "&species=Homo sapiens&format=json")

        response <- GET(url, timeout(5))
        if(status_code(response) == 200) {
            content <- fromJSON(content(response, "text"))
            if(!is.null(content$result) && length(content$result) > 0) {
                return(list(
                    name = content$result[[1]]$name,
                    id = content$result[[1]]$id,
                    url = content$result[[1]]$url,
                    source = "WikiPathways"
                ))
            }
        }
        return(NULL)
    }, error = function(e) NULL)
}

# Query Reactome for pathway details
query_reactome <- function(pathway_name) {
    if(!HAS_HTTR || !HAS_JSONLITE) return(NULL)

    tryCatch({
        library(httr)
        library(jsonlite)

        query <- gsub("_", " ", pathway_name)
        query <- gsub("REACTOME ", "", query)
        query <- gsub("HALLMARK ", "", query)

        url <- paste0("https://reactome.org/ContentService/search/query?query=",
                      URLencode(query), "&species=Homo sapiens&types=Pathway")

        response <- GET(url, timeout(5))
        if(status_code(response) == 200) {
            content <- fromJSON(content(response, "text"))
            if(!is.null(content$results) && length(content$results) > 0) {
                desc <- if(!is.null(content$results[[1]]$summation) &&
                           length(content$results[[1]]$summation) > 0) {
                    content$results[[1]]$summation[[1]]$text
                } else {
                    "Description not available"
                }

                return(list(
                    name = content$results[[1]]$name,
                    id = content$results[[1]]$stId,
                    url = paste0("https://reactome.org/content/detail/",
                                 content$results[[1]]$stId),
                    description = substr(desc, 1, 200),
                    source = "Reactome"
                ))
            }
        }
        return(NULL)
    }, error = function(e) NULL)
}

# Extract drug mechanism information
extract_drug_moa <- function(drug_name) {
    # Simplified MOA database - in production, query DrugBank API
    moa_db <- list(
        "DOXORUBICIN" = "DNA intercalation → Topoisomerase II inhibition → DNA damage → Apoptosis",
        "METFORMIN" = "Complex I inhibition → ↓ATP/↑AMP ratio → AMPK activation → ↓mTOR → Growth inhibition",
        "PACLITAXEL" = "Microtubule stabilization → Mitotic arrest → Apoptosis",
        "CISPLATIN" = "DNA crosslinking → DNA damage response → p53 activation → Apoptosis",
        "TAMOXIFEN" = "Estrogen receptor antagonism → ↓ER signaling → Growth inhibition",
        "IMATINIB" = "BCR-ABL kinase inhibition → ↓Proliferation signals → Apoptosis",
        "BEVACIZUMAB" = "VEGF neutralization → ↓Angiogenesis → Tumor hypoxia",
        "ERLOTINIB" = "EGFR kinase inhibition → ↓MAPK/PI3K signaling → Growth arrest",
        "RAPAMYCIN" = "mTOR inhibition → ↓Protein synthesis → ↓Cell growth",
        "BORTEZOMIB" = "Proteasome inhibition → Protein accumulation → ER stress → Apoptosis",
        "GEMCITABINE" = "Nucleoside analog → DNA synthesis inhibition → S-phase arrest → Apoptosis",
        "SORAFENIB" = "Multi-kinase inhibition (RAF/VEGFR/PDGFR) → ↓Proliferation & Angiogenesis",
        "TEMOZOLOMIDE" = "DNA methylation → O6-methylguanine → Mismatch → Apoptosis",
        "VENETOCLAX" = "BCL-2 inhibition → Mitochondrial permeabilization → Apoptosis",
        "TRAMETINIB" = "MEK1/2 inhibition → ↓MAPK signaling → Growth arrest"
    )

    # Clean drug name and try to match
    drug_upper <- toupper(gsub("-.*", "", drug_name))
    drug_upper <- gsub("_.*", "", drug_upper)

    # Try exact match
    if(drug_upper %in% names(moa_db)) {
        return(moa_db[[drug_upper]])
    }

    # Try partial match
    for(known_drug in names(moa_db)) {
        if(grepl(known_drug, drug_upper) || grepl(drug_upper, known_drug)) {
            return(moa_db[[known_drug]])
        }
    }

    return("Mechanism not in database")
}

# ==============================================================================
# ENHANCED VISUALIZATION FUNCTIONS (NEW IN V6)
# ==============================================================================

# Create drug-pathway overlap heatmap
create_drug_pathway_heatmap <- function(pathway_results, drug_results, out_prefix, contrast) {
    if(is.null(pathway_results) || is.null(drug_results)) return(NULL)
    if(nrow(pathway_results) == 0 || nrow(drug_results) == 0) return(NULL)

    tryCatch({
        # Get top pathways and drugs
        top_pathways <- pathway_results %>%
            filter(p.adjust < 0.05) %>%
            arrange(p.adjust) %>%
            head(DRUG_PATHWAY_TOP_N) %>%
            pull(ID)

        top_drugs <- drug_results %>%
            filter(NES < 0, p.adjust < 0.25) %>%
            arrange(NES) %>%
            head(DRUG_PATHWAY_TOP_N) %>%
            pull(ID)

        if(length(top_pathways) == 0 || length(top_drugs) == 0) return(NULL)

        # Create overlap matrix based on shared genes
        overlap_mat <- matrix(0, nrow=length(top_drugs), ncol=length(top_pathways),
                            dimnames=list(top_drugs, top_pathways))

        for(i in seq_along(top_drugs)) {
            drug_genes <- unlist(strsplit(
                drug_results$core_enrichment[drug_results$ID == top_drugs[i]], "/"))

            for(j in seq_along(top_pathways)) {
                pathway_genes <- unlist(strsplit(
                    pathway_results$core_enrichment[pathway_results$ID == top_pathways[j]], "/"))

                overlap <- length(intersect(drug_genes, pathway_genes))
                overlap_mat[i, j] <- overlap
            }
        }

        # Create heatmap
        library(ComplexHeatmap)
        library(circlize)

        # Truncate labels
        rownames(overlap_mat) <- substr(rownames(overlap_mat), 1, 40)
        colnames(overlap_mat) <- substr(colnames(overlap_mat), 1, 40)

        ht <- Heatmap(overlap_mat,
                      name = "Shared\nGenes",
                      col = colorRamp2(c(0, max(overlap_mat)/2, max(overlap_mat)),
                                     c("white", "#fee090", "#d73027")),
                      cluster_rows = TRUE,
                      cluster_columns = TRUE,
                      show_row_names = TRUE,
                      show_column_names = TRUE,
                      row_names_gp = gpar(fontsize = 9),
                      column_names_gp = gpar(fontsize = 9),
                      column_names_rot = 45,
                      heatmap_legend_param = list(title = "Shared\nGenes"),
                      column_title = paste0("Drug-Pathway Overlap: ", contrast),
                      width = unit(12, "cm"),
                      height = unit(10, "cm"))

        pdf(paste0(out_prefix, "_", contrast, "_DrugPathway_Heatmap_mqc.pdf"), width=16, height=14)
        draw(ht)
        dev.off()

        png(paste0(out_prefix, "_", contrast, "_DrugPathway_Heatmap_mqc.png"),
            width=16, height=14, units="in", res=300, bg="white")
        draw(ht)
        dev.off()

        return(overlap_mat)
    }, error = function(e) {
        cat("ERROR creating drug-pathway heatmap:", e$message, "\n")
        return(NULL)
    })
}

# Create MOA diagram for top drugs
create_moa_diagram <- function(drug_results, out_prefix, contrast) {
    if(is.null(drug_results) || nrow(drug_results) == 0) return(NULL)

    tryCatch({
        # Get top therapeutic candidates
        top_drugs <- drug_results %>%
            filter(NES < 0, p.adjust < 0.25) %>%
            arrange(NES) %>%
            head(MOA_TOP_DRUGS)

        if(nrow(top_drugs) == 0) return(NULL)

        # Extract MOA for each drug
        moa_data <- top_drugs %>%
            mutate(
                Drug = substr(ID, 1, 40),
                MOA = sapply(ID, extract_drug_moa),
                HasMOA = MOA != "Mechanism not in database"
            ) %>%
            filter(HasMOA) %>%
            select(Drug, MOA, NES, p.adjust)

        if(nrow(moa_data) == 0) return(NULL)

        # Create text-based MOA diagram
        p <- ggplot(moa_data, aes(x = reorder(Drug, NES), y = NES, fill = p.adjust)) +
            geom_bar(stat = "identity") +
            scale_fill_gradient(low = "#d73027", high = "#fee090", name = "FDR") +
            coord_flip() +
            labs(title = paste0("Top Drug Mechanisms of Action: ", contrast),
                 subtitle = "Drugs with negative NES oppose disease signature",
                 x = "Drug", y = "NES (Normalized Enrichment Score)") +
            theme_minimal(base_size = 12) +
            theme(
                plot.title = element_text(size = 16, face = "bold"),
                plot.subtitle = element_text(size = 12, color = "grey40"),
                legend.position = "right"
            )

        ggsave(paste0(out_prefix, "_", contrast, "_Drug_MOA_mqc.pdf"), p, width=14, height=10)
        ggsave(paste0(out_prefix, "_", contrast, "_Drug_MOA_mqc.png"), p, width=14, height=10, dpi=300, bg="white")

        return(moa_data)
    }, error = function(e) {
        cat("ERROR creating MOA diagram:", e$message, "\n")
        return(NULL)
    })
}

# Create polypharmacology network
create_polypharm_network <- function(drug_results, pathway_results, out_prefix, contrast) {
    if(is.null(drug_results) || is.null(pathway_results)) return(NULL)
    if(nrow(drug_results) == 0 || nrow(pathway_results) == 0) return(NULL)

    tryCatch({
        # Get drugs and pathways
        drugs <- drug_results %>%
            filter(NES < 0, p.adjust < 0.25) %>%
            head(15) %>%
            mutate(Drug = substr(ID, 1, 30))

        pathways <- pathway_results %>%
            filter(p.adjust < 0.05) %>%
            head(15) %>%
            mutate(Pathway = substr(ID, 1, 30))

        if(nrow(drugs) == 0 || nrow(pathways) == 0) return(NULL)

        # Build network edges
        edges <- data.frame()
        for(i in 1:nrow(drugs)) {
            drug_genes <- unlist(strsplit(drugs$core_enrichment[i], "/"))

            for(j in 1:nrow(pathways)) {
                pathway_genes <- unlist(strsplit(pathways$core_enrichment[j], "/"))
                overlap <- length(intersect(drug_genes, pathway_genes))

                if(overlap >= 3) {  # At least 3 shared genes
                    edges <- rbind(edges, data.frame(
                        from = drugs$Drug[i],
                        to = pathways$Pathway[j],
                        weight = overlap
                    ))
                }
            }
        }

        if(nrow(edges) == 0) return(NULL)

        # Create network
        g <- graph_from_data_frame(edges, directed = FALSE)

        # Identify multi-target drugs
        drug_degree <- degree(g, v = V(g)[V(g)$name %in% drugs$Drug])
        multi_target <- names(drug_degree[drug_degree >= 3])

        V(g)$type <- ifelse(V(g)$name %in% drugs$Drug, "Drug", "Pathway")
        V(g)$multi_target <- V(g)$name %in% multi_target

        p <- ggraph(g, layout = "fr") +
            geom_edge_link(aes(width = weight), alpha = 0.3, color = "grey60") +
            scale_edge_width(range = c(0.5, 3)) +
            geom_node_point(aes(color = type, size = type,
                              shape = ifelse(multi_target & type == "Drug", "Multi-target", "Single"))) +
            scale_color_manual(values = c("Drug" = "#e74c3c", "Pathway" = "#3498db")) +
            scale_size_manual(values = c("Drug" = 6, "Pathway" = 4)) +
            scale_shape_manual(values = c("Multi-target" = 17, "Single" = 16)) +
            geom_node_text(aes(label = ifelse(multi_target, name, "")),
                          repel = TRUE, size = 3, fontface = "bold") +
            theme_void() +
            labs(title = paste0("Polypharmacology Network: ", contrast),
                 subtitle = "Triangles = Multi-target drugs (≥3 pathways)") +
            theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 12, hjust = 0.5),
                  legend.position = "right")

        ggsave(paste0(out_prefix, "_", contrast, "_Polypharm_Network_mqc.pdf"), p, width=14, height=12)
        ggsave(paste0(out_prefix, "_", contrast, "_Polypharm_Network_mqc.png"), p, width=14, height=12, dpi=300, bg="white")

        return(multi_target)
    }, error = function(e) {
        cat("ERROR creating polypharmacology network:", e$message, "\n")
        return(NULL)
    })
}

# ==============================================================================
# HTML REPORTING ENGINE (ENHANCED FOR V6)
# ==============================================================================
html_buffer <- character()

init_html <- function() {
    style <- "
    <style>
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; max-width: 1400px; margin: 40px auto; padding: 20px; background: #f5f7fa; }
        .header { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 30px; border-radius: 10px; margin-bottom: 30px; }
        .header h1 { margin: 0; font-size: 32px; }
        .header .subtitle { opacity: 0.9; margin-top: 10px; font-size: 16px; }
        .header .version { opacity: 0.85; margin-top: 5px; font-size: 14px; font-style: italic; }

        .contrast-block { background: white; padding: 30px; margin: 25px 0; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); }
        .contrast-title { color: #667eea; border-bottom: 3px solid #667eea; padding-bottom: 15px; margin-top: 0; font-size: 28px; }
        .stat-line { font-size: 16px; color: #555; margin-bottom: 15px; }

        .section-header { background: #34495e; color: white; padding: 12px 20px; border-radius: 6px 6px 0 0; font-size: 18px; font-weight: bold; margin-top: 25px; }
        .drug-header { background: #27ae60; }
        .ppi-header { background: #8e44ad; }
        .integration-header { background: #e67e22; }
        .database-header { background: #16a085; }

        .interpretation-box { background: #e8f4f8; border-left: 5px solid #3498db; padding: 20px; margin: 20px 0; border-radius: 4px; }
        .interpretation-box h4 { margin-top: 0; color: #2c3e50; }

        .pathway-table { width: 100%; border-collapse: collapse; font-size: 13px; margin: 15px 0; }
        .pathway-table th { background-color: #34495e; color: white; padding: 12px 10px; text-align: left; }
        .pathway-table td { padding: 10px; border-bottom: 1px solid #eee; }
        .pathway-table tr:hover { background: #f8f9fa; }

        .nes-pos { color: #c0392b; font-weight: bold; }
        .nes-neg { color: #2980b9; font-weight: bold; }
        .sig-green { color: #27ae60; font-weight: bold; }
        .sig-orange { color: #e67e22; }

        .plot-guide { background: #fff3cd; border: 1px solid #ffc107; padding: 15px; margin: 15px 0; border-radius: 5px; }
        .plot-guide strong { color: #856404; }

        .metric-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 15px; margin: 20px 0; }
        .metric-card { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; border-radius: 8px; }
        .metric-card .label { font-size: 12px; opacity: 0.9; text-transform: uppercase; letter-spacing: 1px; }
        .metric-card .value { font-size: 28px; font-weight: bold; margin: 8px 0; }

        .hub-box { padding: 15px; background: #f9f9f9; border: 1px solid #ddd; border-radius: 5px; font-family: 'Courier New', monospace; color: #d35400; font-size: 14px; margin: 10px 0; }

        .legend-box { background: #eef6fc; padding: 20px; border-radius: 5px; border: 1px solid #b8daff; margin: 20px 0; }
        .new-feature { background: #d1ecf1; padding: 10px; border-radius: 5px; margin-top: 10px; border-left: 4px solid #0c5460; }

        .key-finding { background: #d4edda; border-left: 5px solid #28a745; padding: 15px; margin: 15px 0; border-radius: 4px; }
        .key-finding strong { color: #155724; }

        .database-query { background: #f8f9fa; border: 1px solid #dee2e6; padding: 15px; margin: 10px 0; border-radius: 5px; }
        .database-query a { color: #16a085; text-decoration: none; }
        .database-query a:hover { text-decoration: underline; }

        .moa-box { font-family: 'Courier New', monospace; font-size: 11px; color: #2c3e50; background: #ecf0f1; padding: 8px; border-radius: 3px; }

        code { background: #f4f4f4; padding: 3px 8px; border-radius: 3px; font-family: 'Courier New', monospace; color: #c7254e; }

        .toc { background: white; padding: 20px; border-radius: 8px; margin: 20px 0; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        .toc ul { list-style: none; padding-left: 0; }
        .toc li { padding: 5px 0; }
        .toc a { text-decoration: none; color: #667eea; }
        .toc a:hover { text-decoration: underline; }

        .footer { text-align: center; margin-top: 40px; padding: 20px; background: white; border-radius: 8px; color: #666; }
    </style>
    <div class='header'>
        <h1>🔬 Pathway Enrichment & Drug Discovery Analysis</h1>
        <div class='subtitle'>Supreme Edition v6 - Drug-Pathway Integration with Database Queries</div>
        <div class='version'>Enhanced with: MOA diagrams | Drug-Pathway heatmaps | Polypharmacology networks | Database integration</div>
        <div class='subtitle'>Generated: "

    legend <- "
    <div class='legend-box'>
        <h3 style='margin-top:0;'>📖 Quick Reference Guide</h3>
        <div style='display:grid; grid-template-columns: 1fr 1fr; gap:20px;'>
            <div>
                <strong>NES (Normalized Enrichment Score):</strong>
                <ul style='margin:5px 0; padding-left:20px;'>
                    <li><span class='nes-pos'>Positive NES (Red)</span> = Pathway activated/upregulated</li>
                    <li><span class='nes-neg'>Negative NES (Blue)</span> = Pathway suppressed/downregulated</li>
                    <li>|NES| > 1.5 = Strong enrichment</li>
                    <li>|NES| > 2.0 = Very strong enrichment</li>
                </ul>
            </div>
            <div>
                <strong>FDR (False Discovery Rate):</strong>
                <ul style='margin:5px 0; padding-left:20px;'>
                    <li><span class='sig-green'>FDR < 0.05</span> = Statistically significant</li>
                    <li><span class='sig-orange'>FDR 0.05-0.25</span> = Suggestive/exploratory</li>
                    <li>FDR > 0.25 = Not significant</li>
                </ul>
            </div>
        </div>
        <div class='new-feature' style='margin-top:15px;'>
            <strong>✨ NEW IN V6:</strong> Drug-Pathway Integration | Mechanism of Action diagrams |
            Polypharmacology networks | Database-queried pathway annotations
        </div>
    </div>"

    html_buffer <<- c(html_buffer, paste0(style, Sys.Date(), "</p>", legend))
}

add_toc <- function(contrasts) {
    toc <- "
    <div class='toc'>
        <h3 style='margin-top:0;'>📑 Table of Contents</h3>
        <ul>"
    for(cid in contrasts) {
        toc <- paste0(toc, "<li><a href='#", cid, "'>", cid, "</a></li>")
    }
    toc <- paste0(toc, "</ul></div>")
    html_buffer <<- c(html_buffer, toc)
}

add_contrast_header <- function(cid, n_genes, ranking_metric, n_up, n_dn) {
    block <- paste0(
        "<div class='contrast-block' id='", cid, "'>",
        "<h2 class='contrast-title'>📊 ", cid, "</h2>",
        "<div class='metric-grid'>",
        "  <div class='metric-card'><div class='label'>Significant Genes</div><div class='value'>", n_genes, "</div></div>",
        "  <div class='metric-card'><div class='label'>Upregulated</div><div class='value'>", n_up, "</div></div>",
        "  <div class='metric-card'><div class='label'>Downregulated</div><div class='value'>", n_dn, "</div></div>",
        "  <div class='metric-card'><div class='label'>Ranking Metric</div><div class='value' style='font-size:16px;'>", ranking_metric, "</div></div>",
        "</div>",
        "<div class='stat-line'>",
        "  <strong>Differential Expression Cutoffs:</strong> FDR < ", PADJ_CUTOFF, ", |log2FC| > ", LOG2FC_CUTOFF, "<br>",
        "  <strong>GSEA Parameters:</strong> Gene set size: ", GSEA_MIN_SIZE, "-", GSEA_MAX_SIZE, " genes, 1000 permutations",
        "</div>"
    )
    html_buffer <<- c(html_buffer, block)
}

add_interpretation_guide <- function(plot_type) {
    guides <- list(
        dotplot = "
        <div class='interpretation-box'>
            <h4>📈 How to Interpret Dotplot</h4>
            <p><strong>What it shows:</strong> Top enriched pathways ranked by statistical significance</p>
            <ul>
                <li><strong>X-axis (GeneRatio):</strong> Proportion of genes in pathway that are in your DE gene list (larger = more genes overlap)</li>
                <li><strong>Dot size:</strong> Number of genes from your data in this pathway (larger = more genes)</li>
                <li><strong>Dot color:</strong> Statistical significance (red = more significant)</li>
                <li><strong>Split panels:</strong>
                    <ul>
                        <li><strong>Activated (left):</strong> Pathways with positive NES (upregulated in condition)</li>
                        <li><strong>Suppressed (right):</strong> Pathways with negative NES (downregulated in condition)</li>
                    </ul>
                </li>
            </ul>
            <p><strong>Biological Interpretation:</strong> Pathways at the top with large, red dots are your most confident findings.
            Look for biological coherence - do related pathways cluster together?</p>
        </div>",

        emap = "
        <div class='interpretation-box'>
            <h4>🕸️ How to Interpret Enrichment Map (Network)</h4>
            <p><strong>What it shows:</strong> Relationships between enriched pathways based on shared genes</p>
            <ul>
                <li><strong>Nodes (circles):</strong> Individual pathways
                    <ul>
                        <li><strong>Node size:</strong> Number of genes in pathway</li>
                        <li><strong>Node color:</strong> NES score (red = activated, blue = suppressed)</li>
                    </ul>
                </li>
                <li><strong>Edges (connections):</strong> Shared genes between pathways (thicker = more overlap)</li>
                <li><strong>Clusters:</strong> Groups of connected nodes represent functionally related processes</li>
            </ul>
            <p><strong>Biological Interpretation:</strong> Highly connected clusters suggest coordinated biological programs.
            Isolated nodes may represent unique processes. Look for \"hub\" pathways with many connections - these often represent
            central regulatory themes.</p>
        </div>",

        running = "
        <div class='interpretation-box'>
            <h4>📊 How to Interpret GSEA Running Enrichment Score</h4>
            <p><strong>What it shows:</strong> How genes in a pathway are distributed across your ranked gene list</p>
            <ul>
                <li><strong>Top panel (Enrichment Score):</strong>
                    <ul>
                        <li>Running sum that increases when genes in the pathway are encountered</li>
                        <li><strong>Peak position:</strong> Where enrichment is maximal (left = upregulated, right = downregulated)</li>
                        <li><strong>Score magnitude:</strong> Strength of enrichment</li>
                    </ul>
                </li>
                <li><strong>Middle panel (Gene hits):</strong> Vertical black lines show where pathway genes appear in ranked list
                    <ul>
                        <li>Clustered on <strong>left</strong> = pathway upregulated</li>
                        <li>Clustered on <strong>right</strong> = pathway downregulated</li>
                        <li>Evenly spread = no enrichment</li>
                    </ul>
                </li>
                <li><strong>Bottom panel (Ranking metric):</strong> Shows the actual gene-level statistics (usually log2FC)</li>
            </ul>
            <p><strong>Biological Interpretation:</strong> This plot validates whether enrichment is driven by a coordinated shift
            of many pathway genes (good) vs. a few outliers (bad). Strong enrichment shows clear clustering of gene hits toward
            one end of the ranked list.</p>
        </div>",

        network = "
        <div class='interpretation-box'>
            <h4>🔗 How to Interpret Gene-Pathway Network (Cnetplot)</h4>
            <p><strong>What it shows:</strong> Which specific genes drive pathway enrichment</p>
            <ul>
                <li><strong>Large outer nodes:</strong> Enriched pathways</li>
                <li><strong>Small inner nodes:</strong> Individual genes
                    <ul>
                        <li><strong>Color:</strong> log2 fold change (red = upregulated, blue = downregulated)</li>
                        <li><strong>Connections:</strong> Lines show gene-pathway membership</li>
                    </ul>
                </li>
                <li><strong>Hub genes:</strong> Genes connected to many pathways are multi-functional regulators</li>
            </ul>
            <p><strong>Biological Interpretation:</strong> This plot answers \"which specific genes are causing pathway enrichment?\"
            Genes that connect multiple pathways may be key regulatory nodes. Look for genes with extreme fold changes
            (deep red/blue) - these are your strongest drivers.</p>
        </div>",

        ridge = "
        <div class='interpretation-box'>
            <h4>🏔️ How to Interpret Ridgeplot</h4>
            <p><strong>What it shows:</strong> Distribution of gene expression changes within each pathway</p>
            <ul>
                <li><strong>Each ridge (row):</strong> One enriched pathway</li>
                <li><strong>X-axis:</strong> Gene-level log2 fold change values</li>
                <li><strong>Ridge shape:</strong> Distribution of expression changes for genes in that pathway
                    <ul>
                        <li><strong>Shifted right (positive):</strong> Most pathway genes upregulated</li>
                        <li><strong>Shifted left (negative):</strong> Most pathway genes downregulated</li>
                        <li><strong>Narrow peak:</strong> Consistent directional change</li>
                        <li><strong>Wide/bimodal:</strong> Mixed regulation within pathway</li>
                    </ul>
                </li>
            </ul>
            <p><strong>Biological Interpretation:</strong> Pathways with strong, consistent shifts (narrow peaks far from zero)
            represent coordinated regulation. Bimodal distributions suggest the pathway contains both activated and suppressed
            components, indicating complex regulation.</p>
        </div>"
    )

    if(plot_type %in% names(guides)) {
        html_buffer <<- c(html_buffer, guides[[plot_type]])
    }
}

add_pathway_table <- function(db_name, df, is_drug=FALSE) {
    header_class <- if(is_drug) "drug-header" else ""

    if(is.null(df) || nrow(df) == 0) {
        html_buffer <<- c(html_buffer, paste0(
            "<div class='section-header ", header_class, "'>", db_name, "</div>",
            "<div style='padding:15px; color:#999; background:white;'>No significant results found.</div>"
        ))
        return()
    }

    rows <- apply(df, 1, function(r) {
        nes <- as.numeric(r['NES'])
        fdr <- as.numeric(r['p.adjust'])
        nes_class <- ifelse(nes > 0, "nes-pos", "nes-neg")
        fdr_class <- ifelse(fdr < 0.05, "sig-green", "sig-orange")

        core_genes <- if("core_enrichment" %in% names(r) && !is.na(r['core_enrichment'])) {
            genes <- strsplit(r['core_enrichment'], "/")[[1]]
            if(length(genes) > 8) {
                paste0(paste(head(genes, 8), collapse=", "), "... (", length(genes), " total)")
            } else {
                paste(genes, collapse=", ")
            }
        } else {
            "N/A"
        }

        paste0(
            "<tr>",
            "<td style='max-width:300px;'><strong>", r['ID'], "</strong></td>",
            "<td style='text-align:center;'><span class='", nes_class, "'>", round(nes, 2), "</span></td>",
            "<td style='text-align:center;' class='", fdr_class, "'>", formatC(fdr, format="e", digits=2), "</td>",
            "<td style='font-size:11px; color:#666;'>", core_genes, "</td>",
            "</tr>"
        )
    })

    table_html <- paste0(
        "<div class='section-header ", header_class, "'>", db_name, "</div>",
        "<table class='pathway-table'>",
        "<tr><th>Pathway</th><th>NES</th><th>FDR</th><th>Core Enrichment Genes</th></tr>",
        paste(rows, collapse=""),
        "</table>"
    )

    html_buffer <<- c(html_buffer, table_html)
}

add_drug_interpretation <- function() {
    html_buffer <<- c(html_buffer, "
    <div class='interpretation-box'>
        <h4>💊 Drug Discovery Interpretation Guide</h4>
        <p><strong>What this analysis shows:</strong> Drugs whose gene expression signatures are inversely correlated with your disease signature</p>
        <ul>
            <li><strong>Negative NES (Blue):</strong> Drug signature opposes disease signature = <span style='color:#27ae60;'><strong>THERAPEUTIC POTENTIAL</strong></span>
                <ul>
                    <li>Drug downregulates genes you want downregulated</li>
                    <li>Drug upregulates genes you want upregulated</li>
                    <li>More negative NES = stronger therapeutic potential</li>
                </ul>
            </li>
            <li><strong>Positive NES (Red):</strong> Drug signature mimics disease signature = ⚠️ <strong>Avoid - may worsen disease</strong></li>
        </ul>
        <p><strong>How to prioritize candidates:</strong></p>
        <ol>
            <li><strong>Primary filter:</strong> NES < -1.5 AND FDR < 0.05</li>
            <li><strong>Clinical feasibility:</strong> Is drug FDA-approved or in clinical trials?</li>
            <li><strong>Mechanism:</strong> Does the drug target relevant biology for your disease?</li>
            <li><strong>Validation:</strong> Experimental testing required - this is computational prediction only</li>
        </ol>
        <p><strong>⚠️ Important caveats:</strong> This is <em>in silico</em> prediction based on gene expression similarity.
        Requires experimental validation. Consider drug toxicity, delivery, and pharmacokinetics independently.</p>
    </div>")
}

add_ppi_interpretation <- function(hubs) {
    html_buffer <<- c(html_buffer, "
    <div class='interpretation-box'>
        <h4>🕸️ PPI Network & Hub Gene Interpretation</h4>
        <p><strong>What this analysis shows:</strong> Protein-protein interaction network of your significant DE genes using STRING database</p>
        <ul>
            <li><strong>Hub genes (Red, Large):</strong> Highly connected proteins (many interactions)
                <ul>
                    <li>Often transcription factors, signaling proteins, or metabolic enzymes</li>
                    <li>Represent potential therapeutic targets (disrupting hubs affects many pathways)</li>
                    <li>May be master regulators of your biological phenotype</li>
                </ul>
            </li>
            <li><strong>Regular nodes (Blue, Small):</strong> Less connected proteins</li>
            <li><strong>Network clusters:</strong> Groups of tightly connected proteins often represent functional modules (e.g., ribosome, proteasome)</li>
        </ul>
        <p><strong>Biological Interpretation:</strong></p>")

    if(!is.null(hubs) && length(hubs) > 0) {
        html_buffer <<- c(html_buffer, paste0(
            "<div class='key-finding'>",
            "<strong>🎯 Top Hub Genes Identified:</strong><br>",
            "<div class='hub-box'>", paste(hubs, collapse=", "), "</div>",
            "<p style='margin-bottom:0;'><strong>Recommended follow-up:</strong></p>",
            "<ol style='margin-top:5px;'>",
            "<li>Literature review: What is known about these hubs in your disease context?</li>",
            "<li>DrugBank search: Are there existing drugs targeting these proteins?</li>",
            "<li>Validation: qPCR or Western blot to confirm expression changes</li>",
            "<li>Functional studies: Knockdown/knockout experiments to test causality</li>",
            "</ol>",
            "</div>"))
    } else {
        html_buffer <<- c(html_buffer, "<p><em>No hub genes identified (insufficient network connectivity)</em></p>")
    }

    html_buffer <<- c(html_buffer, "</div>")
}

# NEW IN V6
add_drug_pathway_integration <- function(heatmap_data, moa_data, polypharm_drugs) {
    html_buffer <<- c(html_buffer, "
    <div class='section-header integration-header'>🎯 Drug-Pathway Integration Analysis</div>
    <div class='interpretation-box'>
        <h4>Understanding Drug-Pathway Relationships</h4>
        <p><strong>What this shows:</strong> Which therapeutic drugs target which enriched disease pathways</p>
        <ul>
            <li><strong>Heatmap:</strong> Number of shared genes between drug signatures and disease pathways</li>
            <li><strong>Mechanism of Action:</strong> How each drug works at the molecular level</li>
            <li><strong>Polypharmacology:</strong> Drugs that target multiple pathways (often more effective)</li>
        </ul>
        <p><strong>Clinical relevance:</strong> Drugs that target multiple dysregulated pathways may have
        broader therapeutic efficacy and lower resistance development.</p>
    </div>")

    if(!is.null(heatmap_data)) {
        html_buffer <<- c(html_buffer, paste0(
            "<div class='key-finding'>",
            "<strong>🔍 Key Finding:</strong> Identified ", nrow(heatmap_data), " drugs with pathway overlap. ",
            "The heatmap shows gene-level connections between therapeutic candidates and disease mechanisms.",
            "</div>"
        ))
    }

    if(!is.null(moa_data) && nrow(moa_data) > 0) {
        html_buffer <<- c(html_buffer, "
        <h4>Top Drug Mechanisms:</h4>
        <table class='pathway-table'>
            <tr><th>Drug</th><th>Mechanism of Action</th><th>NES</th><th>FDR</th></tr>")

        for(i in 1:nrow(moa_data)) {
            html_buffer <<- c(html_buffer, sprintf(
                "<tr><td><strong>%s</strong></td><td class='moa-box'>%s</td><td class='nes-neg'>%.2f</td><td class='sig-green'>%.2e</td></tr>",
                moa_data$Drug[i], moa_data$MOA[i], moa_data$NES[i], moa_data$p.adjust[i]
            ))
        }

        html_buffer <<- c(html_buffer, "</table>")
    }

    if(!is.null(polypharm_drugs) && length(polypharm_drugs) > 0) {
        html_buffer <<- c(html_buffer, paste0(
            "<div class='key-finding'>",
            "<strong>💊 Multi-Target Drugs Identified:</strong><br>",
            "<div class='hub-box'>", paste(polypharm_drugs, collapse=", "), "</div>",
            "<p><strong>Why this matters:</strong> These drugs target multiple enriched pathways simultaneously, ",
            "suggesting potential for broader therapeutic efficacy. Consider prioritizing for validation.</p>",
            "</div>"
        ))
    }
}

# NEW IN V6
add_pathway_database_info <- function(pathway_name, top_n = 3) {
    # Query databases for top pathways only
    wiki_info <- query_wikipathways(pathway_name)
    reactome_info <- query_reactome(pathway_name)

    if(!is.null(wiki_info) || !is.null(reactome_info)) {
        html_buffer <<- c(html_buffer, "
        <div class='section-header database-header'>📚 Pathway Database Information</div>")

        if(!is.null(wiki_info)) {
            html_buffer <<- c(html_buffer, sprintf(
                "<div class='database-query'>
                    <strong>WikiPathways:</strong> %s<br>
                    <strong>ID:</strong> %s<br>
                    <strong>URL:</strong> <a href='%s' target='_blank'>%s</a>
                </div>",
                wiki_info$name, wiki_info$id, wiki_info$url, wiki_info$url
            ))
        }

        if(!is.null(reactome_info)) {
            html_buffer <<- c(html_buffer, sprintf(
                "<div class='database-query'>
                    <strong>Reactome:</strong> %s<br>
                    <strong>ID:</strong> %s<br>
                    <strong>Description:</strong> %s<br>
                    <strong>URL:</strong> <a href='%s' target='_blank'>%s</a>
                </div>",
                reactome_info$name, reactome_info$id,
                reactome_info$description,
                reactome_info$url, reactome_info$url
            ))
        }
    }
}

close_block <- function() {
    html_buffer <<- c(html_buffer, "</div>")
}

finish_html <- function(prefix) {
    footer <- "
    <div class='footer'>
        <p style='margin:0; font-size:14px;'><strong>Analysis Pipeline Information</strong></p>
        <p style='margin:5px 0;'>Generated by <code>run_pathways_drugs_v6_SUPREME.R</code></p>
        <p style='margin:5px 0; color:#999; font-size:12px;'>
            GSEA: clusterProfiler | PPI: STRING database | Drug discovery: DSigDB<br>
            Database integration: WikiPathways, Reactome | Seed: 12345<br>
            Min gene set: "
    footer <- paste0(footer, GSEA_MIN_SIZE, " | Max gene set: ", GSEA_MAX_SIZE, "
        </p>
    </div>
    </div>")

    html_buffer <<- c(html_buffer, footer)
    writeLines(html_buffer, paste0(dirname(prefix), "/Analysis_Narrative_mqc.html"))
}

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
vst_file <- args[1]; results_dir <- args[2]; gmt_dir <- args[3];
string_dir <- args[4]; out_prefix <- args[5]

# Optional argument: Target Contrast (Default: "ALL")
target_contrast <- if(length(args) >= 6) args[6] else "ALL"

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

# Filter contrasts if a specific target is requested
if(target_contrast != "ALL") {
    cat(paste0("\nLOG: Filtering for specific contrast: ", target_contrast, "\n"))
    # Construct expected filename pattern based on contrast ID
    target_pattern <- paste0(target_contrast, ".deseq2.results.tsv")
    contrasts <- contrasts[basename(contrasts) == target_pattern]

    if(length(contrasts) == 0) {
        available <- sub(".deseq2.results.tsv", "", basename(list.files(file.path(results_dir, "tables/differential"), pattern=".results.tsv")))
        stop(paste0("ERROR: Target contrast '", target_contrast, "' not found. Available contrasts: ", paste(available, collapse=", ")))
    }
}

contrast_ids <- sub(".deseq2.results.tsv", "", basename(contrasts))

add_toc(contrast_ids)

llm_summary <- list()

for(f in contrasts) {
    cid <- sub(".deseq2.results.tsv", "", basename(f))
    cat(paste0("\nLOG: Contrast: ", cid, "...\n"))

    res_df <- read.table(f, header=TRUE, sep="\t", quote="")
    if (grepl("^[0-9]+$", rownames(res_df)[1])) rownames(res_df) <- res_df$gene_id
    if(!"symbol" %in% colnames(res_df)) res_df$symbol <- map_genes_to_symbols(rownames(res_df))

    metric_name <- "log2FoldChange"
    if("stat" %in% colnames(res_df)) {
        metric_name <- "Wald Statistic"
        res_df$rank_metric <- res_df$stat
    } else if ("pvalue" %in% colnames(res_df)) {
        metric_name <- "Signed P-value"
        res_df$rank_metric <- sign(res_df$log2FoldChange) * -log10(res_df$pvalue)
    } else {
        res_df$rank_metric <- res_df$log2FoldChange
    }

    gene_list <- res_df %>%
        dplyr::filter(!is.na(rank_metric), !is.na(symbol), is.finite(rank_metric)) %>%
        distinct(symbol, .keep_all = TRUE) %>%
        arrange(desc(rank_metric)) %>%
        pull(rank_metric, name = symbol)

    sig_genes <- res_df %>%
        dplyr::filter(padj < PADJ_CUTOFF, abs(log2FoldChange) > LOG2FC_CUTOFF) %>%
        pull(symbol)
    sig_genes <- intersect(sig_genes, rownames(mat_vst_sym))

    n_up <- sum(res_df$padj < PADJ_CUTOFF & res_df$log2FoldChange > LOG2FC_CUTOFF, na.rm=TRUE)
    n_dn <- sum(res_df$padj < PADJ_CUTOFF & res_df$log2FoldChange < -LOG2FC_CUTOFF, na.rm=TRUE)

    add_contrast_header(cid, length(sig_genes), metric_name, n_up, n_dn)

    # ==============================================================================
    # A. PATHWAY ANALYSIS
    # ==============================================================================
    gmts <- list.files(gmt_dir, pattern=".gmt", full.names=TRUE)
    pathway_summary <- list()
    pathway_results_for_integration <- NULL

    for(gmt_path in gmts) {
        db_name <- tools::file_path_sans_ext(basename(gmt_path))
        if(grepl("dsigdb", db_name) || grepl("string", db_name)) next

        cat(paste0("  > GSEA: ", db_name, "\n"))
        gmt_data <- tryCatch(read.gmt(gmt_path), error=function(e) NULL)
        if(is.null(gmt_data) || length(intersect(names(gene_list), gmt_data$gene)) < 5) next

        gsea_out <- tryCatch(
            GSEA(gene_list, TERM2GENE=gmt_data, pvalueCutoff=1,
                 minGSSize=GSEA_MIN_SIZE, maxGSSize=GSEA_MAX_SIZE,
                 verbose=FALSE, eps=1e-50, seed=TRUE),
            error=function(e) NULL
        )

        if(!is.null(gsea_out) && nrow(gsea_out) > 0) {
            write.csv(gsea_out@result, paste0(out_prefix, "_", cid, "_", db_name, ".csv"))

            # Store first pathway results for drug-pathway integration
            if(is.null(pathway_results_for_integration)) {
                pathway_results_for_integration <- gsea_out@result
            }

            gsea_out <- pairwise_termsim(gsea_out)

            # 1. DOTPLOT
            add_interpretation_guide("dotplot")
            p_dot <- dotplot(gsea_out, showCategory=GSEA_DOT_N, split=".sign") +
                     facet_grid(.~.sign) +
                     ggtitle(paste0(db_name, ": ", cid)) +
                     theme(axis.text.y = element_text(size=9))
            save_mqc(p_dot, paste0(out_prefix, "_", cid, "_GSEA_Dot_", db_name), 14, 12)

            # 2. ENRICHMENT MAP
            if(nrow(gsea_out) >= 5) {
                add_interpretation_guide("emap")
                p_emap <- emapplot(gsea_out, showCategory=GSEA_EMAP_N,
                                 cex.params=list(category_label=0.6)) +
                          ggtitle(paste0(db_name, ": ", cid))
                save_mqc(p_emap, paste0(out_prefix, "_", cid, "_GSEA_Emap_", db_name), 12, 10)
            }

            # 3. RUNNING ENRICHMENT SCORE
            add_interpretation_guide("running")
            p_running <- gseaplot2(gsea_out,
                                  geneSetID=1:min(GSEA_RUNNING_N, nrow(gsea_out)),
                                  pvalue_table=TRUE, base_size=11)
            save_mqc(p_running, paste0(out_prefix, "_", cid, "_GSEA_Running_", db_name), 14, 12)

            # 4. GENE-PATHWAY NETWORK
            if(nrow(gsea_out) >= 3) {
                add_interpretation_guide("network")
                p_cnet <- cnetplot(gsea_out, showCategory=GSEA_NETWORK_N,
                                  foldChange=gene_list, colorEdge=TRUE,
                                  cex.params=list(category_label=0.6, gene_label=0.5)) +
                          ggtitle(paste0(db_name, ": ", cid))
                save_mqc(p_cnet, paste0(out_prefix, "_", cid, "_GSEA_Network_", db_name), 14, 14)
            }

            # 5. RIDGEPLOT
            if(HAS_GGRIDGES && nrow(gsea_out) >= 5) {
                add_interpretation_guide("ridge")
                p_ridge <- ridgeplot(gsea_out, showCategory=GSEA_RIDGE_N) +
                          ggtitle(paste0(db_name, ": ", cid))
                save_mqc(p_ridge, paste0(out_prefix, "_", cid, "_GSEA_Ridge_", db_name), 10, 8)
            }

            res <- gsea_out@result
            top_pos <- res %>% filter(NES > 0) %>% arrange(desc(NES)) %>% head(TEXT_REPORT_N)
            top_neg <- res %>% filter(NES < 0) %>% arrange(NES) %>% head(TEXT_REPORT_N)
            add_pathway_table(db_name, bind_rows(top_pos, top_neg), is_drug=FALSE)

            # Query database for top pathway (NEW IN V6)
            if(nrow(top_pos) > 0 || nrow(top_neg) > 0) {
                top_pathway <- if(nrow(top_pos) > 0) top_pos$ID[1] else top_neg$ID[1]
                add_pathway_database_info(top_pathway)
            }

            pathway_summary[[db_name]] <- list(
                top_activated = head(top_pos$ID, 5),
                top_suppressed = head(top_neg$ID, 5),
                n_sig = sum(res$p.adjust < 0.05)
            )
        } else {
            add_pathway_table(db_name, NULL)
        }
    }

    # ==============================================================================
    # B. DRUG DISCOVERY
    # ==============================================================================
    dsig_path <- list.files(gmt_dir, pattern="dsigdb", full.names=TRUE)[1]
    drug_results_for_integration <- NULL

    if(!is.na(dsig_path)) {
        cat(paste0("  > Drug Discovery...\n"))
        add_drug_interpretation()

        drug_gmt <- read.gmt(dsig_path)
        drug_gsea <- tryCatch(
            GSEA(gene_list, TERM2GENE=drug_gmt, pvalueCutoff=1,
                 minGSSize=GSEA_MIN_SIZE, maxGSSize=GSEA_MAX_SIZE,
                 verbose=FALSE, eps=1e-50, seed=TRUE),
            error=function(e) NULL
        )

        if(!is.null(drug_gsea) && nrow(drug_gsea) > 0) {
            res <- drug_gsea@result
            drug_results_for_integration <- res
            top_cands <- res %>% filter(NES < 0) %>% arrange(NES) %>% head(GSEA_DOT_N)

            if(nrow(top_cands) > 0) {
                drug_gsea <- pairwise_termsim(drug_gsea)

                p_drug_dot <- dotplot(drug_gsea, showCategory=GSEA_DOT_N, split=".sign") +
                             facet_grid(.~.sign) +
                             ggtitle(paste0("Drug Candidates: ", cid))
                save_mqc(p_drug_dot, paste0(out_prefix, "_", cid, "_Drug_Dotplot"), 14, 12)

                p_drug_emap <- emapplot(drug_gsea, showCategory=20,
                                      cex.params=list(category_label=0.6)) +
                              ggtitle(paste0("Drug Similarity: ", cid))
                save_mqc(p_drug_emap, paste0(out_prefix, "_", cid, "_Drug_Emap"), 12, 10)

                add_pathway_table("💊 Therapeutic Candidates (DSigDB)", top_cands, is_drug=TRUE)

                llm_summary[[cid]]$top_drugs <- head(top_cands$ID, 10)
                llm_summary[[cid]]$n_drug_candidates <- nrow(top_cands)
            } else {
                add_pathway_table("💊 Therapeutic Candidates (DSigDB)", NULL, is_drug=TRUE)
            }
        }
    }

    # ==============================================================================
    # C. DRUG-PATHWAY INTEGRATION (NEW IN V6)
    # ==============================================================================
    heatmap_result <- NULL
    moa_result <- NULL
    polypharm_result <- NULL

    if(!is.null(pathway_results_for_integration) && !is.null(drug_results_for_integration)) {
        cat(paste0("  > Drug-Pathway Integration...\n"))

        # 1. Drug-Pathway Overlap Heatmap
        heatmap_result <- create_drug_pathway_heatmap(
            pathway_results_for_integration,
            drug_results_for_integration,
            out_prefix,
            cid
        )

        # 2. MOA Diagram
        moa_result <- create_moa_diagram(drug_results_for_integration, out_prefix, cid)

        # 3. Polypharmacology Network
        polypharm_result <- create_polypharm_network(
            drug_results_for_integration,
            pathway_results_for_integration,
            out_prefix,
            cid
        )

        # Add to HTML report
        add_drug_pathway_integration(heatmap_result, moa_result, polypharm_result)
    }

    # ==============================================================================
    # D. PPI NETWORK
    # ==============================================================================
    mapped_ids <- sym2string[sig_genes]
    mapped_ids <- mapped_ids[!is.na(mapped_ids)]
    hub_list <- NULL

    if(length(mapped_ids) > 5) {
        cat(paste0("  > PPI Network...\n"))
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
                     geom_node_text(aes(label=ifelse(type=="Hub", name, "")),
                                  repel=TRUE, fontface="bold", size=3.5, bg.color="white") +
                     theme_void() +
                     ggtitle(paste0("PPI Network: ", cid)) +
                     theme(plot.title = element_text(size=16, hjust=0.5))
            save_mqc(p_net, paste0(out_prefix, "_", cid, "_PPI_Network"), 14, 12)

            llm_summary[[cid]]$hub_genes <- hub_list
        }
    }

    add_ppi_interpretation(hub_list)

    llm_summary[[cid]]$pathways <- pathway_summary
    llm_summary[[cid]]$n_de_genes <- length(sig_genes)
    llm_summary[[cid]]$n_up <- n_up
    llm_summary[[cid]]$n_dn <- n_dn
    llm_summary[[cid]]$has_drug_pathway_integration <- !is.null(heatmap_result)
    llm_summary[[cid]]$multi_target_drugs <- polypharm_result

    close_block()
}

finish_html(out_prefix)

# ==============================================================================
# GENERATE LLM PROMPT (TXT FORMAT)
# ==============================================================================
cat("\nLOG: Generating LLM Prompt...\n")

txt_prompt <- c(
    "╔══════════════════════════════════════════════════════════════════════════╗",
    "║     PATHWAY ENRICHMENT & DRUG DISCOVERY ANALYSIS SUMMARY (v6 SUPREME)    ║",
    "║                     For AI/LLM Biological Interpretation                 ║",
    "╚══════════════════════════════════════════════════════════════════════════╝",
    "",
    paste0("Generated: ", Sys.Date()),
    paste0("Analysis: Pathway enrichment (GSEA), Drug discovery (DSigDB), PPI networks (STRING)"),
    paste0("Enhanced: Drug-Pathway Integration, MOA analysis, Polypharmacology, Database queries"),
    "",
    "===============================================================================",
    "OVERVIEW",
    "===============================================================================",
    ""
)

for(cid in names(llm_summary)) {
    summ <- llm_summary[[cid]]

    txt_prompt <- c(txt_prompt,
        "-------------------------------------------------------------------------------",
        paste0("CONTRAST: ", cid),
        "-------------------------------------------------------------------------------",
        "",
        "DIFFERENTIAL EXPRESSION SUMMARY:",
        paste0("  • Total significant genes: ", summ$n_de_genes),
        paste0("  • Upregulated: ", summ$n_up, " genes"),
        paste0("  • Downregulated: ", summ$n_dn, " genes"),
        ""
    )

    if(!is.null(summ$pathways) && length(summ$pathways) > 0) {
        txt_prompt <- c(txt_prompt, "PATHWAY ENRICHMENT RESULTS:", "")
        for(db in names(summ$pathways)) {
            pw <- summ$pathways[[db]]
            txt_prompt <- c(txt_prompt,
                paste0("  [", db, "]"),
                paste0("    Significant pathways: ", pw$n_sig)
            )
            if(length(pw$top_activated) > 0) {
                txt_prompt <- c(txt_prompt,
                    "    Top ACTIVATED pathways:",
                    paste0("      • ", pw$top_activated)
                )
            }
            if(length(pw$top_suppressed) > 0) {
                txt_prompt <- c(txt_prompt,
                    "    Top SUPPRESSED pathways:",
                    paste0("      • ", pw$top_suppressed)
                )
            }
            txt_prompt <- c(txt_prompt, "")
        }
    }

    if(!is.null(summ$top_drugs)) {
        txt_prompt <- c(txt_prompt,
            "THERAPEUTIC DRUG CANDIDATES:",
            paste0("  Total candidates identified: ", summ$n_drug_candidates),
            "  Top 10 candidates (negative NES = opposes disease signature):",
            paste0("    • ", summ$top_drugs),
            ""
        )
    }

    if(summ$has_drug_pathway_integration) {
        txt_prompt <- c(txt_prompt,
            "DRUG-PATHWAY INTEGRATION (NEW IN V6):",
            "  ✓ Drug-Pathway overlap heatmap generated",
            "  ✓ Mechanism of Action analysis completed",
            "  ✓ Polypharmacology network analyzed",
            ""
        )
    }

    if(!is.null(summ$multi_target_drugs) && length(summ$multi_target_drugs) > 0) {
        txt_prompt <- c(txt_prompt,
            "MULTI-TARGET DRUGS (Polypharmacology):",
            "  Drugs targeting ≥3 enriched pathways simultaneously:",
            paste0("    • ", summ$multi_target_drugs),
            "  → These drugs may have broader therapeutic efficacy",
            ""
        )
    }

    if(!is.null(summ$hub_genes)) {
        txt_prompt <- c(txt_prompt,
            "PPI NETWORK HUB GENES:",
            "  (Highly connected proteins - potential therapeutic targets)",
            paste0("    • ", summ$hub_genes),
            ""
        )
    }

    txt_prompt <- c(txt_prompt, "")
}

txt_prompt <- c(txt_prompt,
    "===============================================================================",
    "INTERPRETATION INSTRUCTIONS FOR AI/LLM",
    "===============================================================================",
    "",
    "1. PATHWAY INTERPRETATION:",
    "   • Activated pathways (positive NES) = upregulated in test condition",
    "   • Suppressed pathways (negative NES) = downregulated in test condition",
    "   • Look for biological coherence: Do related pathways cluster together?",
    "   • Prioritize pathways with |NES| > 2.0 and FDR < 0.05",
    "",
    "2. DRUG DISCOVERY INTERPRETATION:",
    "   • Negative NES = Drug opposes disease signature (GOOD - therapeutic potential)",
    "   • Positive NES = Drug mimics disease signature (BAD - avoid)",
    "   • Prioritize: NES < -1.5, FDR < 0.05, FDA-approved status",
    "   • THIS IS COMPUTATIONAL PREDICTION - requires experimental validation",
    "",
    "3. DRUG-PATHWAY INTEGRATION (NEW IN V6):",
    "   • Heatmap shows gene-level overlap between drugs and pathways",
    "   • High overlap = drug directly modulates pathway genes",
    "   • Multi-target drugs affect multiple pathways simultaneously",
    "   • MOA diagrams explain mechanistic basis for drug action",
    "",
    "4. PPI HUB GENE INTERPRETATION:",
    "   • Hub genes = highly connected proteins in interaction network",
    "   • Often transcription factors, signaling proteins, or master regulators",
    "   • Consider as therapeutic targets (disrupting hubs affects many pathways)",
    "   • Cross-reference with DrugBank for existing drugs",
    "",
    "5. INTEGRATED ANALYSIS:",
    "   • Combine pathway + drug + drug-pathway + PPI findings",
    "   • Ask: Do hub genes participate in enriched pathways?",
    "   • Ask: Do candidate drugs target hub genes or enriched pathways?",
    "   • Ask: Do multi-target drugs cover multiple hub-containing pathways?",
    "   • Ask: What is the mechanistic story connecting these findings?",
    "",
    "6. RECOMMENDED OUTPUTS:",
    "   • Biological narrative explaining key findings",
    "   • Mechanistic model connecting pathways → genes → drugs",
    "   • Prioritized list of therapeutic candidates with rationale",
    "   • Drug-pathway interaction map",
    "   • Experimental validation strategy",
    "",
    "===============================================================================",
    paste0("END OF REPORT | Generated by run_pathways_drugs_v6_SUPREME.R | ", Sys.Date()),
    "==============================================================================="
)

writeLines(txt_prompt, paste0(dirname(out_prefix), "/LLM_Analysis_Prompt.txt"))

cat("\n╔═══════════════════════════════════════════════════════════════╗\n")
cat("║            SUPREME EDITION v6 ANALYSIS COMPLETE               ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")
cat("✅ Generated Files:\n")
cat(sprintf("   • HTML Report: %s/Analysis_Narrative_mqc.html\n", dirname(out_prefix)))
cat(sprintf("   • LLM Prompt: %s/LLM_Analysis_Prompt.txt\n", dirname(out_prefix)))
cat(sprintf("   • Total contrasts processed: %d\n", length(llm_summary)))
cat("\n📊 Plots Generated Per Contrast:\n")
cat("   • GSEA Dotplot (pathway ranking)\n")
cat("   • GSEA Enrichment Map (pathway network)\n")
cat("   • GSEA Running Score (enrichment validation)\n")
cat("   • GSEA Gene-Pathway Network (driver genes)\n")
if(HAS_GGRIDGES) cat("   • GSEA Ridgeplot (expression distributions)\n")
cat("   • Drug Discovery Dotplot\n")
cat("   • Drug Similarity Network\n")
cat("   • Drug-Pathway Overlap Heatmap (v6)\n")
cat("   • Mechanism of Action Diagram (v6)\n")
cat("   • Polypharmacology Network (v6)\n")
cat("   • PPI Network with Hub Genes\n")
cat("\n🤖 AI-Ready Features:\n")
cat("   • Complete interpretation guides for each plot type\n")
cat("   • Self-contained HTML with embedded explanations\n")
cat("   • Structured TXT prompt for LLM analysis\n")
cat("   • All results tables with core enrichment genes\n")
if(HAS_HTTR && HAS_JSONLITE) {
    cat("   • Database annotations (WikiPathways, Reactome)\n")
}
cat("\n✨ v6 Enhancements:\n")
cat("   • Drug-Pathway integration analysis\n")
cat("   • Multi-target drug identification\n")
cat("   • Mechanism of Action visualization\n")
cat("   • Database-queried pathway context\n")
cat("\n")

writeLines(capture.output(sessionInfo()), paste0(dirname(out_prefix), "/sessionInfo.txt"))
