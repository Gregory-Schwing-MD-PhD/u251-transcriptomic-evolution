#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH -o logs/prep_website_%j.out
#SBATCH -e logs/prep_website_%j.err

set -e

# --- 1. SETUP ---
WORK_DIR="$(pwd)/work"
OUT_DIR="$(pwd)/ANALYSIS/cibersortx/for_upload"
mkdir -p "$OUT_DIR" "logs"

# Inputs
NEFTEL_EXPR="$(pwd)/ANALYSIS/refs/neftel_2019/IDHwtGBM.processed.SS2.logTPM.txt.gz"
NEFTEL_META="$(pwd)/ANALYSIS/refs/neftel_2019/IDHwt.GBM.Metadata.SS2.txt"
TPM_FILE="$(pwd)/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_tpm.tsv"

# Container
IMG_R="${HOME}/singularity_cache/go2432-bioconductor.sif"
if [[ ! -f "$IMG_R" ]]; then
    mkdir -p "${HOME}/singularity_cache"
    singularity pull "$IMG_R" docker://go2432/bioconductor:latest
fi

# --- 2. RUN PREPARATION ---
echo "Generating CIBERSORTx Input Files..."

# BINDING PWD to /data
# SCRIPT LOCATION: /data/prepare_cibersort.R (Since it's in the base dir now)
singularity exec --bind $PWD:/data --pwd /data "$IMG_R" \
    Rscript prepare_cibersort.R \
    "$NEFTEL_EXPR" \
    "$NEFTEL_META" \
    "$TPM_FILE"

# Move the outputs to the upload folder for clarity
mv refsample.txt "$OUT_DIR/"
mv phenotypes.txt "$OUT_DIR/"
mv mixture.txt "$OUT_DIR/"

echo "DONE."
echo "Upload these 3 files from ${OUT_DIR} to cibersortx.stanford.edu:"
echo "1. refsample.txt"
echo "2. phenotypes.txt"
echo "3. mixture.txt"
