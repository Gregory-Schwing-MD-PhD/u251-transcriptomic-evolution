#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --job-name=cib_frac
#SBATCH -o logs/run_01_fractions_%j.out
#SBATCH -e logs/run_01_fractions_%j.err

set -e

# --- 0. CREDENTIALS ---
if [[ -z "$CIBERSORTX_TOKEN" ]]; then echo "ERROR: Token missing"; exit 1; fi

# --- 1. SETUP ---
WORK_DIR="$(pwd)/work"
OUT_DIR="$(pwd)/ANALYSIS/cibersortx"
mkdir -p "$OUT_DIR/plots" "logs"

# Inputs
NEFTEL_EXPR="$(pwd)/ANALYSIS/refs/neftel_2019/IDHwtGBM.processed.SS2.logTPM.txt.gz"
NEFTEL_META="$(pwd)/ANALYSIS/refs/neftel_2019/IDHwt.GBM.Metadata.SS2.txt"
TPM_FILE="$(pwd)/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_tpm.tsv"
METADATA_FILE="$(pwd)/ANALYSIS/metadata_therapy.csv" 

# Pull Container
export CIBERSORTX_IMG="${HOME}/singularity_cache/cibersortx_fractions.sif"
if [[ ! -f "$CIBERSORTX_IMG" ]]; then
    mkdir -p "${HOME}/singularity_cache"
    singularity pull "$CIBERSORTX_IMG" docker://containers.cibersortx.stanford.edu/cibersortx/fractions:latest
fi

# --- 2. RUN FRACTIONS (Nextflow) ---
echo "STEP 1: Generating Signature & Fractions..."
nextflow run cibersortx.nf \
    -profile singularity \
    --token "$CIBERSORTX_TOKEN" \
    --email "grego@wayne.edu" \
    --neftel_expr "$NEFTEL_EXPR" \
    --neftel_meta "$NEFTEL_META" \
    --u251_tpm "$TPM_FILE" \
    --outdir "$OUT_DIR" \
    -w "${WORK_DIR}" \
    -resume \
    -with-singularity "$CIBERSORTX_IMG"

# --- 3. VISUALIZE FRACTIONS ---
echo "STEP 2: Visualizing Cell States..."
IMG_R="${HOME}/singularity_cache/go2432-bioconductor.sif"

# Note: We bind PWD to /data and run the R script directly from there
singularity exec --bind $PWD:/data --pwd /data "$IMG_R" \
    Rscript plot_cibersort_kitchen_sink.R \
    "${OUT_DIR}/CIBERSORTx_Results.txt" \
    "$METADATA_FILE" \
    "${OUT_DIR}/plots/Fractions_Analysis"

echo "DONE. Signature Matrix saved at: ${OUT_DIR}/CIBERSORTx_Sig.txt"
