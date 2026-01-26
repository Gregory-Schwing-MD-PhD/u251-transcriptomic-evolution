#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
neftel_expr_path <- args[1] 
neftel_meta_path <- args[2] 
u251_tpm_path    <- args[3] 

print("LOG: Starting CIBERSORTx Preparation...")

# --- 1. LOAD METADATA & FIX HEADER ---
if (!file.exists(neftel_meta_path)) stop("Metadata file not found!")

# Read metadata (fill=TRUE handles ragged rows if any)
neftel_meta <- fread(neftel_meta_path, fill=TRUE)

# CRITICAL FIX: Remove the "TYPE" row if it exists
if (neftel_meta[[1]][1] == "TYPE") {
    print("Detected Single Cell Portal 'TYPE' row. Removing it.")
    neftel_meta <- neftel_meta[-1, ]
}

# Ensure the column used for linking is called "NAME"
if (!"NAME" %in% colnames(neftel_meta)) {
    print("WARNING: 'NAME' column not found. Using first column as Sample ID.")
    colnames(neftel_meta)[1] <- "NAME"
}

# Filter for malignant cells (MES, AC, OPC, NPC)
# We look at 'CellAssignment'
malignant_cells <- neftel_meta %>%
  filter(CellAssignment %in% c("MES-like", "AC-like", "OPC-like", "NPC-like")) %>%
  pull(NAME)

print(paste("Malignant Cells Selected:", length(malignant_cells)))

if (length(malignant_cells) == 0) {
    stop("CRITICAL ERROR: No cells matched. Check if 'CellAssignment' column exists.")
}

# --- 2. LOAD EXPRESSION ---
print("LOG: Loading Expression Data...")
# Handle GZ files via system command
if (grepl("\\.gz$", neftel_expr_path)) {
    neftel_expr <- fread(cmd = paste("gunzip -c", neftel_expr_path))
} else {
    neftel_expr <- fread(neftel_expr_path)
}

# Select Genes and Cells
# Force first column to be "Gene"
colnames(neftel_expr)[1] <- "GENE"

# Intersect cells (Ensure we only ask for cells that exist in Expression)
valid_cells <- intersect(malignant_cells, colnames(neftel_expr))
print(paste("Cells found in Expression Matrix:", length(valid_cells)))

if (length(valid_cells) == 0) stop("CRITICAL ERROR: Metadata IDs do not match Expression columns.")

# Subset
neftel_clean <- neftel_expr %>% select(GENE, all_of(valid_cells))

# Convert log2(TPM+1) -> Linear TPM
# Formula: 2^x - 1
mat_data <- as.matrix(neftel_clean[,-1])
mat_linear <- (2^mat_data) - 1
final_ref <- cbind(Gene = neftel_clean$GENE, as.data.frame(mat_linear))

# Write Reference
fwrite(final_ref, "refsample.txt", sep="\t", quote=FALSE)

# Write Phenotypes
phenotypes <- neftel_meta %>%
  filter(NAME %in% valid_cells) %>%
  select(NAME, CellAssignment)
  
fwrite(phenotypes, "phenotypes.txt", sep="\t", col.names=FALSE, quote=FALSE)

print("LOG: Reference Prepared Successfully.")

# --- 3. PREPARE MIXTURE ---
print("LOG: Loading U251 Data...")
u251 <- fread(u251_tpm_path)

# Handle different column names for gene symbol
gene_col <- grep("gene_name|symbol|gene_id", colnames(u251), ignore.case=TRUE, value=TRUE)[1]
if (is.na(gene_col)) gene_col <- colnames(u251)[1] # Fallback to col 1

print(paste("Using U251 Gene Column:", gene_col))
colnames(u251)[which(colnames(u251) == gene_col)] <- "gene_name"

u251_final <- u251 %>%
  select(gene_name, where(is.numeric)) %>%
  group_by(gene_name) %>%
  summarise(across(everything(), sum))

fwrite(u251_final, "mixture.txt", sep="\t", quote=FALSE)
print("LOG: Preparation Complete.")
