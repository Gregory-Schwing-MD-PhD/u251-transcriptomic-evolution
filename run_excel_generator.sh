#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=1:00:00
#SBATCH --job-name=u251_excel
#SBATCH --output=excel_gen_%j.out
#SBATCH --error=excel_gen_%j.err

set -x

# 1. SETUP (Copied exact environment config)
export CONDA_PREFIX="${HOME}/mambaforge/envs/nextflow"
export PATH="${CONDA_PREFIX}/bin:$PATH"
unset JAVA_HOME
nextflow -v

# 2. CACHE
export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
mkdir -p $XDG_RUNTIME_DIR $NXF_SINGULARITY_CACHEDIR

# 3. SAFETY
export NXF_SINGULARITY_HOME_MOUNT=true
unset LD_LIBRARY_PATH
unset PYTHONPATH
unset R_LIBS
unset R_LIBS_USER
unset R_LIBS_SITE

# ==============================================================================
# 4. PATH CONFIGURATION
# ==============================================================================
# Define Inputs
RESULTS_DIR="ANALYSIS/results_evolution"
COUNTS_FILE="${RESULTS_DIR}/tables/abundance/all.counts.tsv"

# The directory where the GSEA CSV files live (from the previous step)
VIZ_DIR="ANALYSIS/results_visualization/Final_Report"

# The Output File
OUT_XLSX="${VIZ_DIR}/U251_Evolution_Master_Data.xlsx"

# Container Image Config
CONTAINER="docker://go2432/bioconductor:latest"
IMG_PATH="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"

# CHECK AND PULL IF MISSING (Fixed logic)
if [[ ! -f "$IMG_PATH" ]]; then
    echo "Container image not found at $IMG_PATH"
    echo "Pulling container from $CONTAINER..."
    singularity pull "$IMG_PATH" "$CONTAINER"
fi

# ==============================================================================
# 5. EXECUTE R EXCEL GENERATOR
# ==============================================================================
echo "RUNNING STEP: GENERATE MASTER EXCEL FILE"

# Verify inputs exist
if [[ ! -f "$COUNTS_FILE" ]]; then
    echo "ERROR: Counts file not found: $COUNTS_FILE"
    # Fallback to alternate location if pipeline moved things
    COUNTS_FILE="${RESULTS_DIR}/tables/counts/all.counts.tsv"
    echo "Trying alternate path: $COUNTS_FILE"
fi

singularity exec --bind $PWD:/data --pwd /data "$IMG_PATH" \
    Rscript generate_excel.R \
    "$COUNTS_FILE" \
    "$RESULTS_DIR" \
    "$VIZ_DIR" \
    "$OUT_XLSX"

echo "--------------------------------------------------------"
echo "EXCEL GENERATION COMPLETE"
echo "File: $OUT_XLSX"
echo "--------------------------------------------------------"
