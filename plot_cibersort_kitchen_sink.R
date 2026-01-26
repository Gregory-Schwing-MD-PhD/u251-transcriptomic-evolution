#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(tidyr); library(readr)
    library(ggpubr); library(RColorBrewer); library(ComplexHeatmap); library(circlize)
    library(tibble)
})

# --- UTILS ---
save_plot <- function(plot_obj, filename_base, w=9, h=8) {
    tryCatch({
        if (inherits(plot_obj, "Heatmap")) {
            pdf(paste0(filename_base, ".pdf"), width=w, height=h)
            draw(plot_obj)
            dev.off()
            png(paste0(filename_base, ".png"), width=w, height=h, units="in", res=300)
            draw(plot_obj)
            dev.off()
        } else {
            ggsave(paste0(filename_base, ".pdf"), plot_obj, width=w, height=h)
            ggsave(paste0(filename_base, ".png"), plot_obj, width=w, height=h, dpi=300, bg="white")
        }
        cat(paste0("SUCCESS: Saved ", basename(filename_base), "\n"))
        return(TRUE)
    }, error = function(e) {
        cat(paste0("ERROR: Failed to save ", basename(filename_base), ": ", e$message, "\n"))
        return(FALSE)
    })
}

# --- MAIN ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: script.R <Results.txt> <Metadata.csv> <OutPrefix>")

cib_path  <- args[1]
meta_path <- args[2]
out_prefix <- args[3]
dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)

neftel_colors <- c("MES-like"="#d73027", "AC-like"="#4575b4", "OPC-like"="#fdae61", "NPC-like"="#abd9e9")

# Load Data
cib_res <- read_tsv(cib_path, show_col_types = FALSE)
colnames(cib_res)[1] <- "Sample"

fraction_cols <- c("MES-like", "AC-like", "OPC-like", "NPC-like")
if (!all(fraction_cols %in% colnames(cib_res))) {
    fraction_cols <- colnames(cib_res)[!colnames(cib_res) %in% c("Sample", "P-value", "Correlation", "RMSE")]
}

if (file.exists(meta_path)) {
    meta <- read_csv(meta_path, show_col_types = FALSE)
    if (!"sample" %in% tolower(colnames(meta))) colnames(meta)[1] <- "sample"
    plot_data <- cib_res %>% inner_join(meta, by = c("Sample" = "sample"))
} else {
    plot_data <- cib_res
    plot_data$Classification <- "Unknown"
}

# Plot 1: Stacked Bar
long_data <- plot_data %>%
    select(Sample, Classification, all_of(fraction_cols)) %>%
    pivot_longer(cols = all_of(fraction_cols), names_to = "CellState", values_to = "Proportion")

p_bar <- ggplot(long_data, aes(x = Sample, y = Proportion, fill = CellState)) +
    geom_bar(stat = "identity", position = "fill", width = 0.8) +
    scale_fill_manual(values = neftel_colors) +
    facet_grid(~Classification, scales = "free_x", space = "free_x") +
    theme_bw(base_size = 14) +
    labs(title = "Cellular State Composition", y = "Fraction", x = "") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
save_plot(p_bar, paste0(out_prefix, "_Composition_Bar_mqc"), 12, 6)

# Plot 2: Boxplots
if (length(unique(plot_data$Classification)) > 1) {
    p_box <- ggplot(long_data, aes(x = Classification, y = Proportion, fill = CellState)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.8) +
        facet_wrap(~CellState, scales = "free_y", ncol = 4) +
        scale_fill_manual(values = neftel_colors) +
        theme_bw(base_size = 14) +
        labs(title = "Shift in Cell States by Condition", y = "Proportion") +
        stat_compare_means(label = "p.signif", method = "t.test", label.y.npc = 0.9)
    save_plot(p_box, paste0(out_prefix, "_Therapy_Boxplot_mqc"), 12, 6)
}

# Plot 3: PCA
pca_mat <- as.matrix(plot_data[, fraction_cols])
rownames(pca_mat) <- plot_data$Sample
pca <- prcomp(pca_mat, scale. = TRUE)
pca_coords <- as.data.frame(pca$x)
pca_coords$Sample <- rownames(pca_coords)
pca_coords <- left_join(pca_coords, plot_data %>% select(Sample, Classification), by = "Sample")
var_expl <- round(summary(pca)$importance[2, 1:2] * 100, 1)

p_pca <- ggplot(pca_coords, aes(x = PC1, y = PC2, fill = Classification)) +
    geom_point(shape = 21, size = 5, color = "black") +
    theme_bw(base_size = 14) +
    labs(title = "PCA of Cell State Fractions", x = paste0("PC1 (", var_expl[1], "%)"), y = paste0("PC2 (", var_expl[2], "%)")) +
    scale_fill_brewer(palette = "Set1")
save_plot(p_pca, paste0(out_prefix, "_CellState_PCA_mqc"), 8, 8)
