#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH -o logs/prep_website_%j.out
#SBATCH -e logs/prep_website_%j.err

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
# 3. SAFETY LOCKS
# ==========================================
export NXF_SINGULARITY_HOME_MOUNT=true
unset LD_LIBRARY_PATH
unset PYTHONPATH
unset R_LIBS
unset R_LIBS_USER
unset R_LIBS_SITE

# ==========================================
# 4. PREPARATION EXECUTION
# ==========================================
OUT_DIR="$(pwd)/ANALYSIS/cibersortx/for_upload"
mkdir -p "$OUT_DIR"

# Inputs
NEFTEL_EXPR="$(pwd)/ANALYSIS/refs/neftel_2019/IDHwtGBM.processed.SS2.logTPM.txt.gz"
NEFTEL_META="$(pwd)/ANALYSIS/refs/neftel_2019/IDHwt.GBM.Metadata.SS2.txt"
TPM_FILE="$(pwd)/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_tpm.tsv"

# Container
IMG_R="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"
if [[ ! -f "$IMG_R" ]]; then
    echo "Pulling container..."
    singularity pull "$IMG_R" docker://go2432/bioconductor:latest
fi

echo "Generating CIBERSORTx Input Files..."

# SCRIPT LOCATION: /data/prepare_cibersort.R (Mapped from PWD)
singularity exec --bind $PWD:/data --pwd /data "$IMG_R" \
    Rscript prepare_cibersort.R \
    "$NEFTEL_EXPR" \
    "$NEFTEL_META" \
    "$TPM_FILE"

# Check if output files exist before moving
if [[ -f "refsample.txt" ]]; then mv refsample.txt "$OUT_DIR/"; fi
if [[ -f "phenotypes.txt" ]]; then mv phenotypes.txt "$OUT_DIR/"; fi
if [[ -f "mixture.txt" ]]; then mv mixture.txt "$OUT_DIR/"; fi

echo "DONE."
echo "Upload these 3 files from ${OUT_DIR} to cibersortx.stanford.edu:"
echo "1. refsample.txt"
echo "2. phenotypes.txt"
echo "3. mixture.txt"
