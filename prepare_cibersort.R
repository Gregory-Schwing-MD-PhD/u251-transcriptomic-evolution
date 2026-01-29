#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
neftel_expr_path <- args[1] 
neftel_meta_path <- args[2] 
u251_tpm_path    <- args[3] 

print("LOG: Starting CIBERSORTx Preparation...")

# --- 1. LOAD METADATA ---
if (!file.exists(neftel_meta_path)) stop("Metadata file not found!")
neftel_meta <- fread(neftel_meta_path, colClasses = "character", fill=TRUE)

# Clean up header/type rows
if (neftel_meta[[1]][1] == "TYPE") {
    neftel_meta <- neftel_meta[-1, ]
}
if (!"NAME" %in% colnames(neftel_meta)) colnames(neftel_meta)[1] <- "NAME"

# *** CRITICAL FIX 1: SANITIZE SAMPLE IDS IN META ***
# Replace hyphens with underscores
neftel_meta$NAME <- str_replace_all(neftel_meta$NAME, "-", "_")

print("LOG: Processing Cell States...")

# Filter for Malignant Cells
tumor_cells <- neftel_meta %>% filter(CellAssignment == "Malignant")
if (nrow(tumor_cells) == 0) stop("CRITICAL: No cells found with CellAssignment == 'Malignant'")

# Score Columns
score_cols <- c("MESlike1", "MESlike2", "AClike", "OPClike", "NPClike1", "NPClike2")
tumor_scores <- tumor_cells %>%
  select(NAME, all_of(score_cols)) %>%
  mutate(across(all_of(score_cols), as.numeric))

# Assign Dominant State
assign_state <- function(r) {
    s_mes <- max(unname(r["MESlike1"]), unname(r["MESlike2"]), na.rm=TRUE)
    s_npc <- max(unname(r["NPClike1"]), unname(r["NPClike2"]), na.rm=TRUE)
    s_ac  <- unname(r["AClike"])
    s_opc <- unname(r["OPClike"])
    
    # *** CRITICAL FIX 2: SANITIZE CLASS NAMES ***
    # Use underscores for classes (AC_like instead of AC-like)
    scores <- c("MES_like"=s_mes, "NPC_like"=s_npc, "AC_like"=s_ac, "OPC_like"=s_opc)
    return(names(which.max(scores)))
}

tumor_cells$ComputedState <- apply(tumor_scores[, -1], 1, assign_state)
valid_cells <- tumor_cells$NAME # These are now underscored

# --- 2. LOAD EXPRESSION & SYNC ---
print("LOG: Loading Expression Data...")
if (grepl("\\.gz$", neftel_expr_path)) {
    neftel_expr <- fread(cmd = paste("gunzip -c", neftel_expr_path))
} else {
    neftel_expr <- fread(neftel_expr_path)
}
colnames(neftel_expr)[1] <- "Gene"

# *** CRITICAL FIX 3: SANITIZE SAMPLE IDS IN EXPRESSION ***
# Force column headers to match the underscored metadata
clean_cols <- str_replace_all(colnames(neftel_expr), "-", "_")
setnames(neftel_expr, colnames(neftel_expr), clean_cols)

# Intersect cells
common_cells <- intersect(valid_cells, colnames(neftel_expr))
print(paste("Matching cells found:", length(common_cells)))

if (length(common_cells) == 0) stop("CRITICAL: No matching cells! Check ID formats.")

# Subset Expression
neftel_clean <- neftel_expr %>% select(Gene, all_of(common_cells))

# Convert to Linear TPM (2^x - 1)
mat_data <- as.matrix(neftel_clean[,-1])
mat_linear <- (2^mat_data) - 1
final_ref <- cbind(Gene = neftel_clean$Gene, as.data.frame(mat_linear))

# Write Reference
fwrite(final_ref, "refsample.txt", sep="\t", quote=FALSE)

# Write Phenotypes
phenotypes <- tumor_cells %>%
  filter(NAME %in% common_cells) %>%
  select(NAME, ComputedState)
  
fwrite(phenotypes, "phenotypes.txt", sep="\t", col.names=FALSE, quote=FALSE)

print("LOG: Reference & Phenotypes Prepared Successfully.")

# --- 3. PREPARE MIXTURE ---
print("LOG: Loading U251 Data...")
u251 <- fread(u251_tpm_path)
cols <- colnames(u251)
if ("gene_name" %in% cols) gene_col <- "gene_name" else gene_col <- cols[1]

u251_final <- u251 %>%
  select(all_of(gene_col), where(is.numeric)) %>%
  rename(Gene = all_of(gene_col)) %>%
  group_by(Gene) %>%
  summarise(across(everything(), sum))

fwrite(u251_final, "mixture.txt", sep="\t", quote=FALSE)
print("LOG: Preparation Complete.")
