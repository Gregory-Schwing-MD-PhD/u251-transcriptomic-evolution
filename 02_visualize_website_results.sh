#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH -o logs/viz_website_%j.out
#SBATCH -e logs/viz_website_%j.err

set -e

# --- 1. SETUP ---
WEBSITE_DIR="$(pwd)/ANALYSIS/cibersortx/website_results"
METADATA_FILE="$(pwd)/ANALYSIS/metadata_therapy.csv"
OUT_DIR="$(pwd)/ANALYSIS/cibersortx/plots"
mkdir -p "$OUT_DIR"

# Container
IMG_R="${HOME}/singularity_cache/go2432-bioconductor.sif"

# --- 2. FIND RESULTS ---
# Auto-detect the Results file (It usually has a random Job ID in the name)
RESULTS_FILE=$(find "$WEBSITE_DIR" -name "*Results.txt" | head -n 1)

if [[ -z "$RESULTS_FILE" ]]; then
    echo "ERROR: Could not find a *Results.txt file in $WEBSITE_DIR"
    echo "Did you unzip the website download?"
    exit 1
fi

echo "Found Results: $RESULTS_FILE"

# --- 3. RUN FRACTIONS VISUALIZATION ---
echo "Visualizing Cell Fractions..."

singularity exec --bind $PWD:/data --pwd /data "$IMG_R" \
    Rscript bin/plot_cibersort_kitchen_sink.R \
    "$RESULTS_FILE" \
    "$METADATA_FILE" \
    "${OUT_DIR}/Website_Fractions"


# --- 4. RUN HIRES VISUALIZATION (If files exist) ---
# Check if there are any HiRes output files (usually start with CIBERSORTxHiRes)
if ls ${WEBSITE_DIR}/CIBERSORTxHiRes*.txt 1> /dev/null 2>&1; then
    echo "Found High-Res files. Running Virtual DE..."
    
    mkdir -p "${OUT_DIR}/hires_volcanoes"
    
    singularity exec --bind $PWD:/data --pwd /data "$IMG_R" \
        Rscript bin/plot_hires_diffexpr.R \
        "$WEBSITE_DIR" \
        "$METADATA_FILE" \
        "${OUT_DIR}/hires_volcanoes"
else
    echo "No High-Res files found. Skipping Step 4."
fi

echo "DONE. Check ${OUT_DIR} for PDFs."
