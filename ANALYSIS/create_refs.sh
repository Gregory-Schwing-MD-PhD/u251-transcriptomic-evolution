#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status.

# ==============================================================================
# SCRIPT: setup_refs.sh (Permission-Proof Version)
# CONTEXT: Run this script INSIDE the 'ANALYSIS' directory.
#          It will create ./refs/human, ./refs/rat, etc.
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

# URLs
RAT_FASTA_URL="http://ftp.ensembl.org/pub/release-110/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz"
RAT_GTF_URL="http://ftp.ensembl.org/pub/release-110/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.110.gtf.gz"

# 1. Process FASTA
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
# SECTION 3: PATHWAY DATABASES
# ==============================================================================
echo -e "\n[3/3] Processing Pathway Databases..."
cd "$PATHWAY_DIR" || exit 1

# 3A. Human Hallmarks (MSigDB)
FILE="h.all.v2023.2.Hs.symbols.gmt"
if [ -f "$FILE" ] && [ -s "$FILE" ]; then
    echo "  [SKIP] $FILE exists."
else
    echo "  [DOWN] Downloading $FILE..."
    wget --no-check-certificate -O "$FILE" "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt"
fi

# 3B. Brain.GMT (GitHub)
# Human
FILE="BrainGMTv2_HumanOrthologs.gmt.txt"
if [ -f "$FILE" ] && [ -s "$FILE" ]; then
    echo "  [SKIP] $FILE exists."
else
    echo "  [DOWN] Downloading $FILE..."
    wget --no-check-certificate -O "$FILE" "https://raw.githubusercontent.com/hagenaue/Brain_GMT/main/BrainGMTv2_HumanOrthologs.gmt.txt"
fi

# Rat
FILE="BrainGMTv2_RatOrthologs.gmt.txt"
if [ -f "$FILE" ] && [ -s "$FILE" ]; then
    echo "  [SKIP] $FILE exists."
else
    echo "  [DOWN] Downloading $FILE..."
    wget --no-check-certificate -O "$FILE" "https://raw.githubusercontent.com/hagenaue/Brain_GMT/main/BrainGMTv2_RatOrthologs.gmt.txt"
fi

# 3C. Generate Rat Hallmarks (R Script)
FILE="r.hallmark.v2023.2.Rn.symbols.gmt"
if [ -f "$FILE" ] && [ -s "$FILE" ]; then
    echo "  [SKIP] Rat Hallmarks ($FILE) already generated."
else
    echo "  [GEN]  Generating Rat Hallmarks via R (msigdbr)..."
    if ! command -v Rscript &> /dev/null; then
        echo "     [WARNING] Rscript not found. Skipping generation."
    else
        # Define local library path for R
        export R_LIBS_USER="$RLIB_DIR"
        
        cat <<EOF > generate_rat_gmt.R
        # Ensure the local library dir exists
        lib_dir <- "$RLIB_DIR"
        if (!dir.exists(lib_dir)) dir.create(lib_dir, recursive = TRUE)
        .libPaths(c(lib_dir, .libPaths()))
        
        # Install msigdbr if missing
        if (!require("msigdbr", quietly = TRUE)) {
            if (!require("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager", repos="http://cran.us.r-project.org", lib=lib_dir)
            }
            BiocManager::install("msigdbr", update=FALSE, ask=FALSE, lib=lib_dir)
        }
        
        library(msigdbr)
        rat_hallmarks <- msigdbr(species = "Rattus norvegicus", category = "H")
        rat_list <- split(x = rat_hallmarks\$gene_symbol, f = rat_hallmarks\$gs_name)
        
        file_conn <- file("$FILE", "w")
        for (pathway in names(rat_list)) {
            genes <- rat_list[[pathway]]
            line <- paste(c(pathway, paste0("https://www.gsea-msigdb.org/gsea/msigdb/cards/", pathway), genes), collapse = "\t")
            writeLines(line, file_conn)
        }
        close(file_conn)
EOF
        Rscript generate_rat_gmt.R
        rm generate_rat_gmt.R
    fi
fi

echo "========================================================"
echo "Setup Complete!"
echo "References located at: $BASE_DIR"
echo "========================================================"
