#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=2:00:00
#SBATCH --job-name=u251_evolution
#SBATCH -o evolution_%j.out
#SBATCH -e evolution_%j.err

set -x

# 1. SETUP
export CONDA_PREFIX="${HOME}/mambaforge/envs/nextflow"
export PATH="${CONDA_PREFIX}/bin:$PATH"
unset JAVA_HOME

# 2. CACHE
export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
mkdir -p $XDG_RUNTIME_DIR $NXF_SINGULARITY_CACHEDIR

# 3. SAFETY
export NXF_SINGULARITY_HOME_MOUNT=true
unset LD_LIBRARY_PATH
unset R_LIBS

# 4. PATHS
CONTAINER="docker://go2432/bioconductor:latest"
IMG_PATH="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"

if [[ ! -f "$IMG_PATH" ]]; then
    singularity pull "$IMG_PATH" "$CONTAINER"
fi

RESULTS_DIR="ANALYSIS/results_evolution"
METADATA_FILE="ANALYSIS/metadata.csv"

# --- OUTPUT CONFIGURATION ---
# New directory to ensure no overwrites of previous reports
OUT_DIR="ANALYSIS/results_global_subtypes"
# This prefix + _mqc.png is what MultiQC looks for
OUTPUT_PREFIX="${OUT_DIR}/Global_Subtypes"

# 5. EXECUTE R ANALYSIS
echo "RUNNING GLOBAL SUBTYPE & EVOLUTION ANALYSIS..."
singularity exec --bind $PWD:/data --pwd /data "$IMG_PATH" \
    Rscript -e "options(warn=1); source('run_global_subtypes.R')" \
    "$RESULTS_DIR" \
    "$OUTPUT_PREFIX" \
    "$METADATA_FILE"

# 6. RUN MULTIQC
echo "RUNNING FINAL MULTIQC..."
MULTIQC_CONTAINER="docker://multiqc/multiqc:v1.33"

# Points explicitly to OUT_DIR where plots are saved
singularity exec --bind $PWD:/data --pwd /data $MULTIQC_CONTAINER multiqc \
    "$OUT_DIR" \
    --force \
    --title "U251 Global Evolutionary Report" \
    --filename "U251_Global_Subtypes_Report_MultiQC.html" \
    --outdir "$OUT_DIR"

echo "DONE."
