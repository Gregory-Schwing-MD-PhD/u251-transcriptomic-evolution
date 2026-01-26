#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH -o logs/viz_website_%j.out
#SBATCH -e logs/viz_website_%j.err

set -x

# ==========================================
# 1. ENVIRONMENT SETUP (CRITICAL FIX)
# ==========================================
echo "LOG: Setting PATH directly..."
export CONDA_PREFIX="${HOME}/mambaforge/envs/nextflow"
export PATH="${CONDA_PREFIX}/bin:$PATH"
unset JAVA_HOME

# ==========================================
# 2. CACHE & STORAGE SETUP
# ==========================================
export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
mkdir -p $XDG_RUNTIME_DIR $NXF_SINGULARITY_CACHEDIR
WORK_DIR="$(pwd)/work"

# ==========================================
# 3. VISUALIZATION EXECUTION
# ==========================================
# Input/Output Directories
WEBSITE_DIR="$(pwd)/ANALYSIS/cibersortx/website_results"
METADATA_FILE="$(pwd)/ANALYSIS/metadata_therapy.csv"
OUT_DIR="$(pwd)/ANALYSIS/cibersortx/plots"
mkdir -p "$OUT_DIR" "logs"

# Container
IMG_R="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"
if [[ ! -f "$IMG_R" ]]; then
    echo "Pulling container..."
    singularity pull "$IMG_R" docker://go2432/bioconductor:latest
fi

# --- STEP A: FIND RESULTS ---
# Auto-detect the Results file (It usually has a random Job ID in the name)
RESULTS_FILE=$(find "$WEBSITE_DIR" -name "*Results.txt" | head -n 1)

if [[ -z "$RESULTS_FILE" ]]; then
    echo "ERROR: Could not find a *Results.txt file in $WEBSITE_DIR"
    echo "Please unzip your CIBERSORTx download and place the text files in:"
    echo "$WEBSITE_DIR"
    exit 1
fi

echo "Found Results: $RESULTS_FILE"

# --- STEP B: VISUALIZE FRACTIONS ---
echo "Visualizing Cell Fractions..."

# BINDING PWD to /data
# SCRIPT LOCATION: /data/plot_cibersort_kitchen_sink.R
singularity exec --bind $PWD:/data --pwd /data "$IMG_R" \
    Rscript plot_cibersort_kitchen_sink.R \
    "$RESULTS_FILE" \
    "$METADATA_FILE" \
    "${OUT_DIR}/Website_Fractions"


# --- STEP C: VISUALIZE HIRES (If files exist) ---
# Check if there are any HiRes output files (usually start with CIBERSORTxHiRes)
if ls ${WEBSITE_DIR}/CIBERSORTxHiRes*.txt 1> /dev/null 2>&1; then
    echo "Found High-Res files. Running Virtual DE..."
    
    mkdir -p "${OUT_DIR}/hires_volcanoes"
    
    singularity exec --bind $PWD:/data --pwd /data "$IMG_R" \
        Rscript plot_hires_diffexpr.R \
        "$WEBSITE_DIR" \
        "$METADATA_FILE" \
        "${OUT_DIR}/hires_volcanoes"
else
    echo "No High-Res files found. Skipping Step C."
fi

echo "DONE. Check ${OUT_DIR} for PDFs."
