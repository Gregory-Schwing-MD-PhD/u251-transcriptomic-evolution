#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status.

# ==============================================================================
# SCRIPT: setup_refs.sh (Final Manuscript Version)
# CONTEXT: Run this script INSIDE the 'ANALYSIS' directory.
#          It creates ./refs/human, ./refs/rat, and merges specific gene sets.
# ==============================================================================

# 1. Capture the Absolute Root (Current Directory = ANALYSIS)
ROOT_DIR=$(pwd)
BASE_DIR="$ROOT_DIR/refs"
HUMAN_DIR="$BASE_DIR/human"
RAT_DIR="$BASE_DIR/rat"
PATHWAY_DIR="$BASE_DIR/pathways"
RLIB_DIR="$BASE_DIR/R_libs"  # Local library for R packages

# 2. Create Directory Structure
mkdir -p "$HUMAN_DIR" "$RAT_DIR" "$PATHWAY_DIR" "$RLIB_DIR"

echo "========================================================"
echo "Starting Reference Setup"
echo "Running inside: $ROOT_DIR"
echo "Target output:  $BASE_DIR"
echo "R Library:      $RLIB_DIR"
echo "========================================================"

# ------------------------------------------------------------------------------
# HELPER FUNCTIONS
# ------------------------------------------------------------------------------

check_md5() {
    local file=$1
    local md5_ref=$2
    if [ ! -f "$md5_ref" ]; then return 1; fi
    grep "$file" "$md5_ref" | md5sum -c - --status
    return $?
}

check_gzip() {
    local file=$1
    gzip -t "$file" 2>/dev/null
    return $?
}

# ==============================================================================
# SECTION 1: HUMAN REFERENCE (GENCODE Release 44)
# ==============================================================================
echo -e "\n[1/3] Processing Human Reference (GRCh38 - GENCODE v44)..."
cd "$HUMAN_DIR" || exit 1

# URLs
HUMAN_FASTA_URL="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz"
HUMAN_GTF_URL="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz"
HUMAN_MD5_URL="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/MD5SUMS"

# 1. Get MD5
wget --no-check-certificate -N "$HUMAN_MD5_URL"

# 2. Process FASTA
FILE="GRCh38.primary_assembly.genome.fa.gz"
if [ -f "$FILE" ] && check_md5 "$FILE" "MD5SUMS"; then
    echo "  [SKIP] $FILE exists and MD5 verified."
else
    echo "  [DOWN] Downloading $FILE..."
    wget --no-check-certificate -O "$FILE" "$HUMAN_FASTA_URL"
    echo "  [VFY]  Verifying MD5..."
    check_md5 "$FILE" "MD5SUMS" || { echo "  [ERR] MD5 mismatch for $FILE"; exit 1; }
fi

# 3. Process GTF
FILE="gencode.v44.primary_assembly.annotation.gtf.gz"
if [ -f "$FILE" ] && check_md5 "$FILE" "MD5SUMS"; then
    echo "  [SKIP] $FILE exists and MD5 verified."
else
    echo "  [DOWN] Downloading $FILE..."
    wget --no-check-certificate -O "$FILE" "$HUMAN_GTF_URL"
    echo "  [VFY]  Verifying MD5..."
    check_md5 "$FILE" "MD5SUMS" || { echo "  [ERR] MD5 mismatch for $FILE"; exit 1; }
fi

# 4. Create symlink for naming consistency
if [ ! -f "GRCh38.primary_assembly.annotation.gtf.gz" ]; then
    ln -s "$FILE" "GRCh38.primary_assembly.annotation.gtf.gz"
    echo "  [LINK] Created symlink for consistency."
fi

# ==============================================================================
# SECTION 2: RAT REFERENCE (Ensembl Release 110)
# ==============================================================================
echo -e "\n[2/3] Processing Rat Reference (mRatBN7.2 - Ensembl v110)..."
cd "$RAT_DIR" || exit 1

# URLs - UPDATED TO TOPLEVEL (Fixes 404 Error)
RAT_FASTA_URL="http://ftp.ensembl.org/pub/release-110/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz"
RAT_GTF_URL="http://ftp.ensembl.org/pub/release-110/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.110.gtf.gz"

# 1. Process FASTA (Updated filename)
FILE="Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz"
if [ -f "$FILE" ] && check_gzip "$FILE"; then
    echo "  [SKIP] $FILE exists and GZIP integrity OK."
else
    echo "  [DOWN] Downloading $FILE..."
    wget --no-check-certificate -O "$FILE" "$RAT_FASTA_URL"
    echo "  [VFY]  Verifying GZIP integrity..."
    check_gzip "$FILE" || { echo "  [ERR] GZIP integrity failed for $FILE"; exit 1; }
fi

# 2. Process GTF
FILE="Rattus_norvegicus.mRatBN7.2.110.gtf.gz"
if [ -f "$FILE" ] && check_gzip "$FILE"; then
    echo "  [SKIP] $FILE exists and GZIP integrity OK."
else
    echo "  [DOWN] Downloading $FILE..."
    wget --no-check-certificate -O "$FILE" "$RAT_GTF_URL"
    echo "  [VFY]  Verifying GZIP integrity..."
    check_gzip "$FILE" || { echo "  [ERR] GZIP integrity failed for $FILE"; exit 1; }
fi

# ==============================================================================
# SECTION 3: PATHWAY DATABASES (Auto-Clean & Merge)
# ==============================================================================
echo -e "\n[3/3] Processing Pathway Databases..."
cd "$PATHWAY_DIR" || exit 1

# --- Define Filenames ---
HUMAN_HALLMARK="h.all.v2023.2.Hs.symbols.gmt"
HUMAN_C2="c2.all.v2023.2.Hs.symbols.gmt"
HUMAN_BRAIN="BrainGMTv2_HumanOrthologs.gmt.txt"
HUMAN_COMBINED="combined_human.gmt"

# Note: This file will contain generated Hallmark AND GO:BP for Rat
RAT_GENERATED="r.hallmark_and_gobp.v2023.2.Rn.symbols.gmt"
RAT_BRAIN="BrainGMTv2_RatOrthologs.gmt.txt"
RAT_COMBINED="combined_rat.gmt"

# 3A. Human Hallmarks (MSigDB)
if [ -f "$HUMAN_HALLMARK" ] && [ -s "$HUMAN_HALLMARK" ]; then
    echo "  [SKIP] Human Hallmark exists."
else
    echo "  [DOWN] Downloading Human Hallmark..."
    wget --no-check-certificate -O "$HUMAN_HALLMARK" "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt"
fi

# 3B. Human C2 Curated (GBM Subtypes + Reactome)
if [ -f "$HUMAN_C2" ] && [ -s "$HUMAN_C2" ]; then
    echo "  [SKIP] Human C2 (Curated) exists."
else
    echo "  [DOWN] Downloading Human C2 (Curated)..."
    wget --no-check-certificate -O "$HUMAN_C2" "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c2.all.v2023.2.Hs.symbols.gmt"
fi

# 3C. Brain.GMT (GitHub) - Human
if [ -f "$HUMAN_BRAIN" ] && [ -s "$HUMAN_BRAIN" ]; then
    echo "  [SKIP] Human BrainGMT exists."
else
    echo "  [DOWN] Downloading Human BrainGMT..."
    wget --no-check-certificate -O "$HUMAN_BRAIN" "https://raw.githubusercontent.com/hagenaue/Brain_GMT/main/BrainGMTv2_HumanOrthologs.gmt.txt"
fi

# 3D. Brain.GMT (GitHub) - Rat
if [ -f "$RAT_BRAIN" ] && [ -s "$RAT_BRAIN" ]; then
    echo "  [SKIP] Rat BrainGMT exists."
else
    echo "  [DOWN] Downloading Rat BrainGMT..."
    wget --no-check-certificate -O "$RAT_BRAIN" "https://raw.githubusercontent.com/hagenaue/Brain_GMT/main/BrainGMTv2_RatOrthologs.gmt.txt"
fi

# 3E. CLEANING STEP (Critical for BrainGMT)
echo "  [CLEAN] Sanitizing BrainGMT files (removing Windows line endings and trailing tabs)..."
sed -i 's/[ \t\r]*$//' "$HUMAN_BRAIN"
sed -i 's/[ \t\r]*$//' "$RAT_BRAIN"

# 3F. Generate Rat Gene Sets (Hallmark + GO:BP) via R
if [ -f "$RAT_GENERATED" ] && [ -s "$RAT_GENERATED" ]; then
    echo "  [SKIP] Rat Generated Sets (Hallmark + GO:BP) already exist."
else
    echo "  [GEN]  Generating Rat Sets via R (msigdbr)..."
    if ! command -v Rscript &> /dev/null; then
        echo "       [WARNING] Rscript not found. Skipping generation."
    else
        export R_LIBS_USER="$RLIB_DIR"
        cat <<INNER > generate_rat_gmt.R
        lib_dir <- "$RLIB_DIR"
        if (!dir.exists(lib_dir)) dir.create(lib_dir, recursive = TRUE)
        .libPaths(c(lib_dir, .libPaths()))

        # Install msigdbr
        if (!require("msigdbr", quietly = TRUE)) {
            if (!require("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager", repos="http://cran.us.r-project.org", lib=lib_dir)
            }
            BiocManager::install("msigdbr", update=FALSE, ask=FALSE, lib=lib_dir)
        }
        library(msigdbr)

        # 1. Hallmark (H) - General Inflammation
        rat_hallmarks <- msigdbr(species = "Rattus norvegicus", category = "H")
        rat_list_H <- split(x = rat_hallmarks\$gene_symbol, f = rat_hallmarks\$gs_name)

        # 2. GO Biological Process (C5:BP) - Specific Scarring/Wound Healing
        rat_gobp <- msigdbr(species = "Rattus norvegicus", category = "C5", subcategory = "GO:BP")
        rat_list_BP <- split(x = rat_gobp\$gene_symbol, f = rat_gobp\$gs_name)

        # 3. Combine
        rat_full_list <- c(rat_list_H, rat_list_BP)

        # 4. Write
        file_conn <- file("$RAT_GENERATED", "w")
        for (pathway in names(rat_full_list)) {
            genes <- rat_full_list[[pathway]]
            line <- paste(c(pathway, "http://www.gsea-msigdb.org/gsea/msigdb/cards/", genes), collapse = "\t")
            writeLines(line, file_conn)
        }
        close(file_conn)
INNER
        Rscript generate_rat_gmt.R
        rm generate_rat_gmt.R
    fi
fi

# 3G. MERGE STEP (With Duplicate Removal)
echo "  [MERGE] Concatenating Gene Sets and removing duplicates..."

# --- Human Merge ---
# Priority: Hallmark > BrainGMT > C2 (awk keeps the first instance it finds)
cat "$HUMAN_HALLMARK" "$HUMAN_BRAIN" "$HUMAN_C2" > "$HUMAN_COMBINED.tmp"
awk '!seen[$1]++' "$HUMAN_COMBINED.tmp" > "$HUMAN_COMBINED"
rm "$HUMAN_COMBINED.tmp"

# Count unique sets for the log
HUMAN_COUNT=$(wc -l < "$HUMAN_COMBINED")
echo "      -> Created $HUMAN_COMBINED ($HUMAN_COUNT unique sets)"

# --- Rat Merge ---
# Priority: Hallmark/GOBP > BrainGMT
if [ -f "$RAT_GENERATED" ]; then
    cat "$RAT_GENERATED" "$RAT_BRAIN" > "$RAT_COMBINED.tmp"
    awk '!seen[$1]++' "$RAT_COMBINED.tmp" > "$RAT_COMBINED"
    rm "$RAT_COMBINED.tmp"

    RAT_COUNT=$(wc -l < "$RAT_COMBINED")
    echo "      -> Created $RAT_COMBINED ($RAT_COUNT unique sets)"
else
    echo "      [WARN] Rat Generated Sets missing. Skipping Rat merge."
fi

echo "========================================================"
echo "Setup Complete!"
echo "Combined Human GMT: $PATHWAY_DIR/$HUMAN_COMBINED"
echo "Combined Rat GMT:   $PATHWAY_DIR/$RAT_COMBINED"
echo "========================================================"
