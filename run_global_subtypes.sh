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
#SBATCH --mail-user=go2432@wayne.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -x

# 1. SETUP
export CONDA_PREFIX="${HOME}/mambaforge/envs/nextflow"
export PATH="${CONDA_PREFIX}/bin:$PATH"
unset JAVA_HOME
nextflow -v

# 2. CACHE
export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
mkdir -p $XDG_RUNTIME_DIR $NXF_SINGULARITY_CACHEDIR
WORK_DIR="$(pwd)/work"

# 3. SAFETY
export NXF_SINGULARITY_HOME_MOUNT=true
unset LD_LIBRARY_PATH
unset PYTHONPATH
unset R_LIBS
unset R_LIBS_USER
unset R_LIBS_SITE

# 5. VISUALIZATION (GLOBAL SUBTYPES)
echo "RUNNING STEP 4: GLOBAL SUBTYPES ANALYSIS"

CONTAINER="docker://go2432/bioconductor:latest"
IMG_PATH="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"

# --- UPDATED PATHS ---
RESULTS_DIR="ANALYSIS/results_evolution"
# CHANGED: New directory to avoid overwriting visualziation plots
OUTPUT_PREFIX="ANALYSIS/results_global_subtypes/Global_Subtypes"
METADATA_FILE="ANALYSIS/metadata.csv"
# (Contrasts and GMT not used by this specific R script, but variables kept for consistency)
CONTRASTS_FILE="ANALYSIS/contrasts.csv"
GMT_FILE="ANALYSIS/refs/pathways/combined_human.gmt"

if [[ ! -f "$IMG_PATH" ]]; then
    singularity pull "$IMG_PATH" "$CONTAINER"
fi

# R Arguments: <results_dir> <out_prefix> <meta>
# CHANGED: Runs 'run_global_subtypes.R' with only the 3 required arguments
singularity exec --bind $PWD:/data --pwd /data "$IMG_PATH" \
    Rscript -e "options(warn=1); source('run_global_subtypes.R')" \
    "$RESULTS_DIR" \
    "$OUTPUT_PREFIX" \
    "$METADATA_FILE"

# 6. MULTIQC
echo "RUNNING STEP 5: FINAL MULTIQC"
cat << 'CONFIG' > mqc_config.yaml
disable_version_detection: true
sections:
  software_versions:
    hide: true
custom_content:
  order:
    - pathway_analysis
CONFIG

MULTIQC_CONTAINER="docker://multiqc/multiqc:v1.33"

# CHANGED: Input dir points to new results; Filename changed to avoid overwrite
singularity exec --bind $PWD:/data --pwd /data $MULTIQC_CONTAINER multiqc \
    /data/ANALYSIS/results_global_subtypes \
    /data/ANALYSIS/results_human_final \
    /data/ANALYSIS/xengsort_out \
    --force \
    --config /data/mqc_config.yaml \
    --title "U251 Global Subtypes Analysis" \
    --filename "U251_Global_Subtypes_Report.html" \
    --outdir "/data/ANALYSIS/results_evolution"

echo "DONE."
