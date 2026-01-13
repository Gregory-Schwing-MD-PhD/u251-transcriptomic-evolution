#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
metadata_file <- args[2]
out_prefix <- args[3]

# 1. Load Data
counts <- fread(counts_file)
meta <- fread(metadata_file)
# Filter strictly for the conditions we defined
meta <- meta[condition %in% c("Primary", "Recurrent", "Control_Brain")]

# Match samples
valid_samples <- intersect(meta$sample, colnames(counts))
counts_mat <- as.matrix(counts[, ..valid_samples])
rownames(counts_mat) <- counts$gene_id
meta <- meta[match(valid_samples, meta$sample),]

# 2. Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(counts_mat), colData = meta, design = ~ condition)
dds <- DESeq(dds)

# 3. Get Contrasts
res_litt <- results(dds, contrast=c("condition", "Recurrent", "Primary"))
res_bg   <- results(dds, contrast=c("condition", "Control_Brain", "Primary"))

# 4. Subtraction Logic
df <- data.frame(
    Gene = rownames(res_litt),
    LFC_LITT = res_litt$log2FoldChange,
    Padj_LITT = res_litt$padj,
    LFC_BG = res_bg$log2FoldChange
)
df <- na.omit(df)

df$Class <- "Noise"
df$Class[df$Padj_LITT < 0.05 & df$LFC_LITT > 1.0] <- "LITT_Candidate"
# The Filter: If Normal Brain also has high expression (LFC_BG > 1), it's an artifact
df$Class[df$Class == "LITT_Candidate" & df$LFC_BG > 1.0] <- "Purity_Artifact"

clean_genes <- df[df$Class == "LITT_Candidate", ]
write.csv(clean_genes, paste0(out_prefix, "_Clean_LITT_Genes.csv"), row.names=FALSE)

pdf(paste0(out_prefix, "_Subtraction_Scatter.pdf"), width=8, height=8)
ggplot(df, aes(x=LFC_BG, y=LFC_LITT, color=Class)) +
    geom_point(alpha=0.6) +
    geom_hline(yintercept=1, linetype="dashed") + geom_vline(xintercept=1, linetype="dashed") +
    labs(x="Log2FC (Control Brain vs Primary)", y="Log2FC (Recurrent vs Primary)", title="LITT Signature Subtraction") +
    scale_color_manual(values=c("LITT_Candidate"="red", "Purity_Artifact"="green", "Noise"="gray")) +
    theme_minimal()
dev.off()
