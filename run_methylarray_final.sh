#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=go2432@wayne.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o methyl_prod_%j.out
#SBATCH -e methyl_prod_%j.err

# Trace every command to the log file for debugging
set -x

# ==========================================
# 1. ENVIRONMENT SETUP
# ==========================================
echo "LOG: Setting PATH directly..."
export CONDA_PREFIX="${HOME}/mambaforge/envs/nextflow"
export PATH="${CONDA_PREFIX}/bin:$PATH"
unset JAVA_HOME

echo "LOG: Verifying Nextflow..."
nextflow -v

# ==========================================
# 2. CACHE & STORAGE SETUP
# ==========================================
export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
mkdir -p $XDG_RUNTIME_DIR $NXF_SINGULARITY_CACHEDIR

# Use a specific work directory for methylation
WORK_DIR="$(pwd)/work_methyl"

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
# 4. PIPELINE EXECUTION
# ==========================================

# --- NF-CORE/METHYLARRAY ---
echo "RUNNING: NF-CORE/METHYLARRAY"

nextflow run nf-core/methylarray \
    -profile singularity \
    -params-file methylarray_params.yaml \
    -w "${WORK_DIR}" \
    -resume \
    -ansi-log false

unset NXF_PARAMS

# --- FINAL REPORTING ---
echo "RUNNING: FINAL MULTIQC AGGREGATION"

singularity exec docker://multiqc/multiqc:v1.33 multiqc \
    --force \
    --title "U251 Methylation Analysis" \
    --filename "U251_Methylation_Report.html" \
    --outdir "ANALYSIS/results_methylarray" \
    --exclude software_versions \
    ANALYSIS/results_methylarray

echo "DONE. Final report: ANALYSIS/results_methylarray/U251_Methylation_Report.html"
