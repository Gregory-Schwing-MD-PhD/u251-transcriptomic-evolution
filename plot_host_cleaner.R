#!/usr/bin/env Rscript

# ==========================================
# HOST PURITY FILTER (SIGNAL SUBTRACTION) - RAT/MICROENV
# ==========================================
# UPDATES:
# - Swapped org.Rn.eg.db for EnsDb.Rnorvegicus.v79 (Native Ensembl Rat DB)

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(DESeq2)
    # library(org.Rn.eg.db) # <-- REMOVED
    library(EnsDb.Rnorvegicus.v79) # <-- ADDED
    library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript plot_host_cleaner.R <counts> <metadata> <out_prefix>")

counts_file <- args[1]
metadata_file <- args[2]
out_prefix <- args[3]

# Ensure output directory exists
dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)

# ==========================================
# 1. LOAD DATA
# ==========================================
counts <- fread(counts_file)
meta <- fread(metadata_file)

# UPDATE: Ensure we are using the correct condition names from your metadata
target_conditions <- c("Primary_Microenv", "Recurrent_Microenv", "Control_Brain")
meta <- meta[condition %in% target_conditions]

# Match samples between counts and metadata
valid_samples <- intersect(meta$sample, colnames(counts))
if(length(valid_samples) == 0) stop("ERROR: No matching samples found between counts and metadata.")

counts_mat <- as.matrix(counts[, ..valid_samples])
rownames(counts_mat) <- counts$gene_id
meta <- meta[match(valid_samples, meta$sample),]

# ==========================================
# 2. RUN DESEQ2
# ==========================================
dds <- DESeqDataSetFromMatrix(countData = round(counts_mat), colData = meta, design = ~ condition)
dds <- DESeq(dds)

# ==========================================
# 3. DEFINE CONTRASTS
# ==========================================
# Contrast 1: The Treatment Effect (Recurrent vs Primary)
res_litt <- results(dds, contrast=c("condition", "Recurrent_Microenv", "Primary_Microenv"))

# Contrast 2: The "Infiltration" Effect (Control Brain vs Primary)
res_bg   <- results(dds, contrast=c("condition", "Control_Brain", "Primary_Microenv"))

# ==========================================
# 4. SUBTRACTION LOGIC & ID MAPPING
# ==========================================
df <- data.frame(
    Gene = rownames(res_litt),
    LFC_LITT = res_litt$log2FoldChange,
    Padj_LITT = res_litt$padj,
    LFC_BG = res_bg$log2FoldChange
)
df <- na.omit(df)

# --- MAP SYMBOLS (ENSEMBL RAT) ---
clean_ids <- sub("\\..*", "", df$Gene)

# Map using EnsDb.Rnorvegicus (Keytype is GENEID)
df$Symbol <- suppressWarnings(mapIds(EnsDb.Rnorvegicus.v79, 
                                     keys=clean_ids, 
                                     column="SYMBOL", 
                                     keytype="GENEID", 
                                     multiVals="first"))

# Fallback to ID if symbol not found
df$Symbol[is.na(df$Symbol)] <- df$Gene[is.na(df$Symbol)]

# --- CLASSIFY ---
df$Class <- "Noise"

# Candidate: Significant change in LITT
df$Class[df$Padj_LITT < 0.05 & abs(df$LFC_LITT) > 1.0] <- "LITT_Candidate"

# Artifact: If the gene follows the "Normal Brain" trend too closely
# Case A: Upregulated in LITT, but also Upregulated in Control Brain -> Likely Normal Brain Infiltration
df$Class[df$Class == "LITT_Candidate" & df$LFC_LITT > 0 & df$LFC_BG > 1.0] <- "Purity_Artifact"

# Case B: Downregulated in LITT, but also Downregulated in Control Brain -> Likely Tumor Purity Artifact
df$Class[df$Class == "LITT_Candidate" & df$LFC_LITT < 0 & df$LFC_BG < -1.0] <- "Purity_Artifact"

# ==========================================
# 5. OUTPUTS
# ==========================================
# Save the Clean List (Candidates Only)
clean_genes <- df[df$Class == "LITT_Candidate", ]
write.csv(clean_genes, paste0(out_prefix, "_Clean_LITT_Genes.csv"), row.names=FALSE)
cat(sprintf("LOG: Identified %d clean LITT-response genes (removed %d artifacts).\n",
            nrow(clean_genes), sum(df$Class == "Purity_Artifact")))

# Plot
pdf(paste0(out_prefix, "_Subtraction_Scatter.pdf"), width=8, height=8)
p <- ggplot(df, aes(x=LFC_BG, y=LFC_LITT, color=Class)) +
    geom_point(alpha=0.6) +
    geom_hline(yintercept=c(-1, 1), linetype="dashed") +
    geom_vline(xintercept=c(-1, 1), linetype="dashed") +
    labs(x="Log2FC (Control Brain vs Primary)",
         y="Log2FC (Recurrent vs Primary)",
         title="Signal Subtraction (Rat Microenvironment)",
         subtitle="Removing Normal Brain Purity Artifacts") +
    scale_color_manual(values=c("LITT_Candidate"="red", "Purity_Artifact"="green", "Noise"="gray")) +
    theme_minimal()

# Label top 10 Artifacts (Genes that tricked us!)
top_artifacts <- head(df[df$Class == "Purity_Artifact" & df$LFC_LITT > 2, ], 10)
p <- p + geom_text_repel(data=top_artifacts, aes(label=Symbol), color="darkgreen", max.overlaps=10)

print(p)
dev.off()
