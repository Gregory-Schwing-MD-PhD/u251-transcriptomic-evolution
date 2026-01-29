#!/bin/bash
set -e

# ==============================================================================
# SCRIPT: create_refs_final.sh (Includes STRING v12.0)
# CONTEXT: Run INSIDE 'ANALYSIS'.
# ==============================================================================

ROOT_DIR=$(pwd)
BASE_DIR="$ROOT_DIR/refs"
HUMAN_DIR="$BASE_DIR/human"
RAT_DIR="$BASE_DIR/rat"
PATHWAY_ROOT="$BASE_DIR/pathways"
PATHWAY_HUMAN="$PATHWAY_ROOT/human"
PATHWAY_RAT="$PATHWAY_ROOT/rat"
STRING_DIR="$PATHWAY_HUMAN/human_string"  # NEW: For PPI Networks
RLIB_DIR="$BASE_DIR/R_libs"

mkdir -p "$HUMAN_DIR" "$RAT_DIR" "$PATHWAY_HUMAN" "$PATHWAY_RAT" "$STRING_DIR" "$RLIB_DIR"

echo "========================================================"
echo "Starting Reference Setup"
echo "Human Pathways:   $PATHWAY_HUMAN"
echo "STRING PPI DB:    $STRING_DIR"
echo "Rat Pathways:     $PATHWAY_RAT"
echo "========================================================"

# --- HELPER FUNCTIONS ---
check_md5() {
    local file=$1
    local md5_ref=$2
    if [ ! -f "$md5_ref" ]; then return 1; fi
    grep "$file" "$md5_ref" | md5sum -c - --status
    return $?
}
check_gzip() {
    gzip -t "$1" 2>/dev/null
    return $?
}
check_valid_file() {
    local file=$1
    if [ ! -s "$file" ]; then
        echo "    [ERR] Download failed or file is empty: $file"
        rm -f "$file"
        return 1
    fi
    return 0
}

# ==============================================================================
# SECTIONS 1 & 2: GENOMES
# ==============================================================================
echo -e "\n[1/4] Processing Human Reference (GRCh38 - GENCODE v44)..."
cd "$HUMAN_DIR" || exit 1
wget --no-check-certificate -N "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/MD5SUMS"

FILE="GRCh38.primary_assembly.genome.fa.gz"
if [ -f "$FILE" ] && check_md5 "$FILE" "MD5SUMS"; then
    echo "  [SKIP] $FILE exists and MD5 verified."
else
    echo "  [DOWN] Downloading $FILE..."
    wget --no-check-certificate -O "$FILE" "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz"
    check_md5 "$FILE" "MD5SUMS" || { echo "  [ERR] MD5 mismatch"; exit 1; }
fi

FILE="gencode.v44.primary_assembly.annotation.gtf.gz"
if [ -f "$FILE" ] && check_md5 "$FILE" "MD5SUMS"; then
    echo "  [SKIP] $FILE exists and MD5 verified."
else
    echo "  [DOWN] Downloading $FILE..."
    wget --no-check-certificate -O "$FILE" "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz"
    check_md5 "$FILE" "MD5SUMS" || { echo "  [ERR] MD5 mismatch"; exit 1; }
fi
[ ! -f "GRCh38.primary_assembly.annotation.gtf.gz" ] && ln -s "$FILE" "GRCh38.primary_assembly.annotation.gtf.gz"

echo -e "\n[2/4] Processing Rat Reference (mRatBN7.2 - Ensembl v110)..."
cd "$RAT_DIR" || exit 1

FILE="Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz"
if [ -f "$FILE" ] && check_gzip "$FILE"; then
    echo "  [SKIP] $FILE exists and GZIP integrity OK."
else
    echo "  [DOWN] Downloading $FILE..."
    wget --no-check-certificate -O "$FILE" "http://ftp.ensembl.org/pub/release-110/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz"
    check_gzip "$FILE" || { echo "  [ERR] GZIP check failed"; exit 1; }
fi

FILE="Rattus_norvegicus.mRatBN7.2.110.gtf.gz"
if [ -f "$FILE" ] && check_gzip "$FILE"; then
    echo "  [SKIP] $FILE exists and GZIP integrity OK."
else
    echo "  [DOWN] Downloading $FILE..."
    wget --no-check-certificate -O "$FILE" "http://ftp.ensembl.org/pub/release-110/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.110.gtf.gz"
    check_gzip "$FILE" || { echo "  [ERR] GZIP check failed"; exit 1; }
fi

# ==============================================================================
# SECTION 3: PATHWAYS
# ==============================================================================
echo -e "\n[3/4] Processing Pathway Databases..."

# --- HUMAN ---
cd "$PATHWAY_HUMAN" || exit 1
H_HALLMARK="hs.hallmark.v2023.2.gmt"
H_C2="hs.c2_curated.v2023.2.gmt"
H_BRAIN="hs.brain_gmt_v2.gmt"
H_DBSIG="hs.dsigdb.v1.0.gmt"
H_GOBP="hs.gobp.v2023.2.gmt"
H_KEGG="hs.kegg.v2023.2.gmt"

if [ ! -s "$H_HALLMARK" ]; then
    wget -O "$H_HALLMARK" "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt"
    check_valid_file "$H_HALLMARK" || exit 1
else echo "  [SKIP] Human Hallmark exists."; fi

if [ ! -s "$H_C2" ]; then
    wget -O "$H_C2" "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c2.all.v2023.2.Hs.symbols.gmt"
    check_valid_file "$H_C2" || exit 1
else echo "  [SKIP] Human C2 exists."; fi

if [ ! -s "$H_BRAIN" ]; then
    wget -O "$H_BRAIN" "https://raw.githubusercontent.com/hagenaue/Brain_GMT/main/BrainGMTv2_HumanOrthologs.gmt.txt"
    check_valid_file "$H_BRAIN" || exit 1
    sed -i 's/[ \t\r]*$//' "$H_BRAIN"
else echo "  [SKIP] Human BrainGMT exists."; fi

if [ ! -s "$H_DBSIG" ]; then
    wget -O "$H_DBSIG" "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=DSigDB"
    check_valid_file "$H_DBSIG" || exit 1
else echo "  [SKIP] Human DSigDB exists."; fi

if [ ! -s "$H_GOBP" ]; then
    wget -O "$H_GOBP" "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c5.go.bp.v2023.2.Hs.symbols.gmt"
    check_valid_file "$H_GOBP" || exit 1
else echo "  [SKIP] Human GO:BP exists."; fi

if [ ! -s "$H_KEGG" ]; then
    [ -s "$H_C2" ] || { echo "    [ERR] C2 missing, cannot extract KEGG"; exit 1; }
    grep "^KEGG_" "$H_C2" > "$H_KEGG"
    check_valid_file "$H_KEGG" || exit 1
else echo "  [SKIP] Human KEGG exists."; fi


# --- RAT ---
cd "$PATHWAY_RAT" || exit 1
R_BRAIN="rn.brain_gmt_v2.gmt"
R_HALLMARK="rn.hallmark.v2023.2.gmt"
R_GOBP="rn.gobp.v2023.2.gmt"
R_KEGG="rn.kegg.v2023.2.gmt"

if [ ! -s "$R_BRAIN" ]; then
    wget -O "$R_BRAIN" "https://raw.githubusercontent.com/hagenaue/Brain_GMT/main/BrainGMTv2_RatOrthologs.gmt.txt"
    check_valid_file "$R_BRAIN" || exit 1
    sed -i 's/[ \t\r]*$//' "$R_BRAIN"
else echo "  [SKIP] Rat BrainGMT exists."; fi

# Generate Rat sets via R (Robust to Version Differences)
if [[ ! -s "$R_HALLMARK" || ! -s "$R_GOBP" || ! -s "$R_KEGG" ]]; then
    echo "  [GEN]  Generating Rat sets via R..."

    export R_LIBS_USER="$RLIB_DIR"

    cat <<INNER > generate_rat_sets.R
    lib_dir <- "$RLIB_DIR"
    if (!dir.exists(lib_dir)) dir.create(lib_dir, recursive = TRUE)
    .libPaths(c(lib_dir, .libPaths()))

    if (!require("msigdbr", quietly = TRUE)) {
         if (!require("BiocManager", quietly = TRUE)) {
             install.packages("BiocManager", repos="http://cran.us.r-project.org", lib=lib_dir)
         }
         BiocManager::install("msigdbr", update=FALSE, ask=FALSE, lib=lib_dir)
    }
    library(msigdbr)

    # 1. Load Data
    all_rat <- msigdbr(species = "Rattus norvegicus")
    cols <- colnames(all_rat)

    # 2. Dynamic Column Mapping
    if ("gs_collection" %in% cols) {
        cat_col <- "gs_collection"
        sub_col <- "gs_subcollection"
    } else if ("gs_cat" %in% cols) {
        cat_col <- "gs_cat"
        sub_col <- "gs_subcat"
    } else {
        stop("CRITICAL: Neither 'gs_collection' nor 'gs_cat' found in msigdbr columns.")
    }

    # 3. Write Function
    write_gmt <- function(df, filename) {
        if (nrow(df) == 0) {
            message("WARNING: No rows found for ", filename)
            return()
        }
        gene_list <- split(x = df\$gene_symbol, f = df\$gs_name)
        file_conn <- file(filename, "w")
        for (pathway in names(gene_list)) {
            genes <- gene_list[[pathway]]
            line <- paste(c(pathway, "http://msigdb_ortholog_generated", genes), collapse = "\t")
            writeLines(line, file_conn)
        }
        close(file_conn)
        message(paste("Generated:", filename))
    }

    # 4. Generate Sets
    write_gmt(all_rat[all_rat[[cat_col]] == "H", ], "$R_HALLMARK")
    write_gmt(all_rat[all_rat[[cat_col]] == "C5" & grepl("GO:BP", all_rat[[sub_col]]), ], "$R_GOBP")
    write_gmt(all_rat[all_rat[[cat_col]] == "C2" & grepl("^KEGG_", all_rat\$gs_name), ], "$R_KEGG")

INNER

    Rscript generate_rat_sets.R || { echo "    [ERR] R script failed"; exit 1; }
    rm generate_rat_sets.R

    check_valid_file "$R_HALLMARK" || exit 1
    check_valid_file "$R_GOBP" || exit 1
    check_valid_file "$R_KEGG" || exit 1
else
    echo "  [SKIP] Rat generated sets exist."
fi

# ==============================================================================
# SECTION 4: STRING PPI NETWORK (v12.0)
# ==============================================================================
echo -e "\n[4/4] Processing STRING Database (v12.0 PPI Network)..."
cd "$STRING_DIR" || exit 1

# 1. Protein Links (The Network)
FILE="9606.protein.links.v12.0.txt.gz"
if [ -f "$FILE" ]; then
    echo "  [SKIP] STRING Links exist."
else
    echo "  [DOWN] Downloading STRING Links (Network Topology)..."
    wget --no-check-certificate -O "$FILE" "https://stringdb-static.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
    check_valid_file "$FILE" || exit 1
fi

# 2. Protein Info (The ID Mapping)
FILE="9606.protein.info.v12.0.txt.gz"
if [ -f "$FILE" ]; then
    echo "  [SKIP] STRING Info exists."
else
    echo "  [DOWN] Downloading STRING Info (ID Mapping)..."
    wget --no-check-certificate -O "$FILE" "https://stringdb-static.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
    check_valid_file "$FILE" || exit 1
fi

echo "========================================================"
echo "Setup Complete!"
echo "STRING DB located at: $STRING_DIR"
echo "========================================================"
