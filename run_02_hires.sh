#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --job-name=cib_hires
#SBATCH -o logs/run_02_hires_%j.out
#SBATCH -e logs/run_02_hires_%j.err

set -e

# --- 0. CREDENTIALS ---
if [[ -z "$CIBERSORTX_TOKEN" ]]; then echo "ERROR: Token missing"; exit 1; fi

# --- 1. SETUP ---
OUT_DIR="$(pwd)/ANALYSIS/cibersortx"
HIRES_DIR="${OUT_DIR}/hires"
mkdir -p "$HIRES_DIR" "logs"

# Inputs
SIG_MATRIX="${OUT_DIR}/CIBERSORTx_Sig.txt"
TPM_FILE="$(pwd)/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_tpm.tsv"
METADATA_FILE="$(pwd)/ANALYSIS/metadata_therapy.csv" 

if [[ ! -f "$SIG_MATRIX" ]]; then
    echo "CRITICAL ERROR: Signature Matrix not found! Run run_01_fractions.sh first."
    exit 1
fi

# Pull Container
export CIBERSORTX_HIRES_IMG="${HOME}/singularity_cache/cibersortx_hires.sif"
if [[ ! -f "$CIBERSORTX_HIRES_IMG" ]]; then
    echo "Pulling HiRes Container..."
    mkdir -p "${HOME}/singularity_cache"
    singularity pull "$CIBERSORTX_HIRES_IMG" docker://containers.cibersortx.stanford.edu/cibersortx/hires:latest
fi

# --- 2. RUN HIGH-RES IMPUTATION ---
CELL_TYPES=("MES-like" "AC-like" "OPC-like" "NPC-like")

echo "Starting High-Resolution Purification..."

for CTYPE in "${CELL_TYPES[@]}"; do
    echo "  > Processing: $CTYPE"
    
    if ls ${HIRES_DIR}/*${CTYPE}*.txt 1> /dev/null 2>&1; then
        echo "    Skipping $CTYPE (Output exists)"
        continue
    fi

    singularity exec --bind $PWD:/src/data --bind $OUT_DIR:/src/outdir "$CIBERSORTX_HIRES_IMG" \
        /src/CIBERSORTxHiRes \
        --username "grego@wayne.edu" \
        --token "$CIBERSORTX_TOKEN" \
        --mixture "$TPM_FILE" \
        --sigmatrix "$SIG_MATRIX" \
        --subset "$CTYPE" \
        --rmbatchBmode TRUE \
        --threads 8 \
        --output_dir /src/outdir/hires
done

# --- 3. DIFFERENTIAL EXPRESSION ---
echo "STEP 3: Running Virtual Differential Expression..."
IMG_R="${HOME}/singularity_cache/go2432-bioconductor.sif"

singularity exec --bind $PWD:/data --pwd /data "$IMG_R" \
    Rscript plot_hires_diffexpr.R \
    "$HIRES_DIR" \
    "$METADATA_FILE" \
    "${HIRES_DIR}/plots"

echo "DONE. Check ${HIRES_DIR}/plots for Volcano plots."
