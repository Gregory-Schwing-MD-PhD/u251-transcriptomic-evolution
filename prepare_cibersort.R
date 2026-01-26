#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(tibble)
    library(stringr) # Added for text cleaning
})

args <- commandArgs(trailingOnly = TRUE)
neftel_expr_path <- args[1] 
neftel_meta_path <- args[2] 
u251_tpm_path    <- args[3] 

print("LOG: Starting CIBERSORTx Preparation...")

# --- 1. LOAD METADATA & CALCULATE STATES ---
if (!file.exists(neftel_meta_path)) stop("Metadata file not found!")

# Read as character to safely handle the "TYPE" row
neftel_meta <- fread(neftel_meta_path, colClasses = "character", fill=TRUE)

# Remove the "TYPE" row if present
if (neftel_meta[[1]][1] == "TYPE") {
    print("Detected 'TYPE' row. Removing it.")
    neftel_meta <- neftel_meta[-1, ]
}

if (!"NAME" %in% colnames(neftel_meta)) colnames(neftel_meta)[1] <- "NAME"

print("LOG: Processing Cell States...")

# Filter for Malignant Cells
tumor_cells <- neftel_meta %>% filter(CellAssignment == "Malignant")
if (nrow(tumor_cells) == 0) stop("CRITICAL: No cells found with CellAssignment == 'Malignant'")

# Score Columns
score_cols <- c("MESlike1", "MESlike2", "AClike", "OPClike", "NPClike1", "NPClike2")
missing_cols <- setdiff(score_cols, colnames(tumor_cells))
if (length(missing_cols) > 0) stop("Cannot calculate states; missing score columns.")

# Convert to numeric
tumor_scores <- tumor_cells %>%
  select(NAME, all_of(score_cols)) %>%
  mutate(across(all_of(score_cols), as.numeric))

# Assign Dominant State
assign_state <- function(r) {
    s_mes <- max(unname(r["MESlike1"]), unname(r["MESlike2"]), na.rm=TRUE)
    s_npc <- max(unname(r["NPClike1"]), unname(r["NPClike2"]), na.rm=TRUE)
    s_ac  <- unname(r["AClike"])
    s_opc <- unname(r["OPClike"])
    
    # CRITICAL FIX: Use underscores, not hyphens
    scores <- c("MES_like"=s_mes, "NPC_like"=s_npc, "AC_like"=s_ac, "OPC_like"=s_opc)
    return(names(which.max(scores)))
}

states <- apply(tumor_scores[, -1], 1, assign_state)
tumor_cells$ComputedState <- states

print("Cell State Distribution:")
print(table(tumor_cells$ComputedState))

valid_cells <- tumor_cells$NAME

# --- 2. LOAD EXPRESSION ---
print("LOG: Loading Expression Data...")
if (grepl("\\.gz$", neftel_expr_path)) {
    neftel_expr <- fread(cmd = paste("gunzip -c", neftel_expr_path))
} else {
    neftel_expr <- fread(neftel_expr_path)
}
colnames(neftel_expr)[1] <- "Gene" # Strict requirement for CIBERSORTx

# Intersect cells
common_cells <- intersect(valid_cells, colnames(neftel_expr))
if (length(common_cells) == 0) stop("CRITICAL: No matching cells between metadata and expression.")

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

print("LOG: Reference Prepared Successfully.")


# --- 3. PREPARE MIXTURE ---
print("LOG: Loading U251 Data...")
u251 <- fread(u251_tpm_path)

# Smart Column Detection
cols <- colnames(u251)
if ("gene_name" %in% cols) {
    gene_col <- "gene_name"
} else if ("symbol" %in% cols) {
    gene_col <- "symbol"
} else if ("gene_id" %in% cols) {
    gene_col <- "gene_id"
} else {
    gene_col <- cols[1]
}

print(paste("Using gene column:", gene_col))

# Robust Selection
u251_final <- u251 %>%
  select(all_of(gene_col), where(is.numeric)) %>%
  rename(Gene = all_of(gene_col)) %>% # Rename to Gene to match Ref
  group_by(Gene) %>%
  summarise(across(everything(), sum))

fwrite(u251_final, "mixture.txt", sep="\t", quote=FALSE)
print("LOG: Preparation Complete.")
