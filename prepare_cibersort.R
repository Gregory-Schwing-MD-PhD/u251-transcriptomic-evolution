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

print("LOG: Loading Neftel Data...")

# --- FIX: Handle .gz files using system command (Bypasses missing R.utils) ---
if (grepl("\\.gz$", neftel_expr_path)) {
    print(paste("Detected GZ file. Using system gunzip for:", neftel_expr_path))
    # 'cmd' argument tells fread to run a shell command instead of opening file directly
    neftel_expr <- fread(cmd = paste("gunzip -c", neftel_expr_path))
} else {
    neftel_expr <- fread(neftel_expr_path)
}

neftel_meta <- fread(neftel_meta_path)

# Filter for malignant cells (MES, AC, OPC, NPC)
malignant_cells <- neftel_meta %>%
  filter(CellAssignment %in% c("MES-like", "AC-like", "OPC-like", "NPC-like")) %>%
  pull(NAME)

# Filter expression matrix
neftel_clean <- neftel_expr %>% select(GENE, all_of(malignant_cells))

# Convert log2(TPM+1) -> Linear TPM
# Formula: 2^x - 1
gene_col <- neftel_clean$GENE
mat_data <- as.matrix(neftel_clean[,-1])
mat_linear <- (2^mat_data) - 1
final_ref <- cbind(Gene = gene_col, as.data.frame(mat_linear))

# Write Reference
fwrite(final_ref, "refsample.txt", sep="\t", quote=FALSE)

# Write Phenotypes (Sample <tab> Class)
phenotypes <- neftel_meta %>%
  filter(NAME %in% malignant_cells) %>%
  select(NAME, CellAssignment)
  
fwrite(phenotypes, "phenotypes.txt", sep="\t", col.names=FALSE, quote=FALSE)

print("LOG: Loading U251 Data...")

# --- 2. PREPARE MIXTURE (Bulk U251) ---
u251 <- fread(u251_tpm_path)

# Aggregate by gene name if needed
if ("gene_name" %in% colnames(u251)) {
    u251_final <- u251 %>%
      select(gene_name, where(is.numeric)) %>%
      group_by(gene_name) %>%
      summarise(across(everything(), sum))
} else {
    print("WARNING: 'gene_name' column not found, using first column as gene symbol.")
    colnames(u251)[1] <- "gene_name"
    u251_final <- u251 %>%
      group_by(gene_name) %>%
      summarise(across(where(is.numeric), sum))
}

# Write Mixture
fwrite(u251_final, "mixture.txt", sep="\t", quote=FALSE)

print("LOG: Preparation Complete.")
