#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(readr); library(tidyr)
    library(limma); library(EnhancedVolcano); library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
hires_dir <- args[1]
meta_path <- args[2]
out_dir   <- args[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

meta <- read_csv(meta_path, show_col_types = FALSE)
if (!"Classification" %in% colnames(meta)) stop("Metadata missing 'Classification' column")

files <- list.files(hires_dir, pattern = "CIBERSORTxHiRes_.*\\.txt", full.names = TRUE)

for (f in files) {
    fname <- basename(f)
    ctype <- NA
    for (t in c("MES-like", "AC-like", "OPC-like", "NPC-like")) {
        if (grepl(t, fname)) ctype <- t
    }
    if (is.na(ctype)) next 

    message(paste("Analyzing Virtual:", ctype))
    expr <- read_tsv(f, show_col_types = FALSE) 
    expr <- expr %>% distinct(GeneSymbol, .keep_all = TRUE) %>% 
        column_to_rownames(var = "GeneSymbol") %>%
        as.matrix()

    valid_samples <- intersect(colnames(expr), meta$sample)
    if (length(valid_samples) < 3) next
    
    expr_clean <- expr[, valid_samples]
    meta_clean <- meta %>% filter(sample %in% valid_samples) %>% arrange(match(sample, valid_samples))
    
    condition <- factor(meta_clean$Classification)
    design <- model.matrix(~0 + condition)
    colnames(design) <- levels(condition)
    
    if (ncol(design) < 2) next

    fit <- lmFit(expr_clean, design)
    
    # Check for specific contrast: Recurrent_U2 vs Primary_U2
    if ("Recurrent_U2" %in% colnames(design) && "Primary_U2" %in% colnames(design)) {
        contr <- makeContrasts(Recurrent_U2 - Primary_U2, levels = design)
        fit2 <- contrasts.fit(fit, contr)
        fit2 <- eBayes(fit2)
        res <- topTable(fit2, number = Inf)
        
        p <- EnhancedVolcano(res, lab = rownames(res), x = 'logFC', y = 'P.Value',
            title = paste(ctype, ": Recurrent vs Primary"),
            subtitle = "Virtual Cell Sorting (CIBERSORTx HiRes)",
            pCutoff = 0.05, FCcutoff = 1.0, pointSize = 2.0, labSize = 3.0)
        ggsave(file.path(out_dir, paste0("Volcano_", ctype, ".pdf")), p, width = 8, height = 8)
        write.csv(res, file.path(out_dir, paste0("DE_Results_", ctype, ".csv")))
    }
}
