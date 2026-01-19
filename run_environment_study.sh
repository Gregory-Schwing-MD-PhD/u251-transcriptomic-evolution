#!/bin/bash
#SBATCH -q secondary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=2:00:00
#SBATCH --job-name=u251_env
#SBATCH -o environment_%j.out
#SBATCH -e environment_%j.err
#SBATCH --mail-user=go2432@wayne.edu
#SBATCH --mail-type=BEGIN,END,FAIL

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
# 4. PIPELINE EXECUTION
# ==========================================

# --- FIX METADATA FORMATTING ---
sed -i 's/\r//g' ANALYSIS/metadata.csv

# --- STEP 3: DIFFERENTIAL ABUNDANCE (ENVIRONMENT) ---
echo "RUNNING STEP 3: DIFFERENTIAL ABUNDANCE (ENVIRONMENT)"

# Using the ENVIRONMENT params file
nextflow run nf-core/differentialabundance \
    -r 1.5.0 \
    -profile singularity \
    -params-file environment_params.yaml \
    -w "${WORK_DIR}" \
    -ansi-log false 

# --- STEP 4: CUSTOM VISUALIZATION ---
echo "RUNNING STEP 4: KITCHEN SINK + ENVIRONMENT"

CONTAINER="docker://go2432/bioconductor:latest"
IMG_PATH="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"

# Updated Paths for Environment Analysis
RESULTS_DIR="ANALYSIS/results_environment"
OUTPUT_PREFIX="ANALYSIS/results_environment/plots/Environment"
COUNTS_FILE="ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv"
METADATA_FILE="ANALYSIS/metadata.csv"
CONTRASTS_FILE="ANALYSIS/contrasts_env.csv" # <--- Critical: This triggers the R script to use Environment mode
GMT_FILE="ANALYSIS/refs/pathways/combined_human.gmt"

if [[ ! -f "$IMG_PATH" ]]; then
    singularity pull "$IMG_PATH" "$CONTAINER"
fi

singularity exec --bind $PWD:/data --pwd /data "$IMG_PATH" Rscript plot_evolution_kitchen_sink.R \
    "$RESULTS_DIR" \
    "$OUTPUT_PREFIX" \
    "$COUNTS_FILE" \
    "$METADATA_FILE" \
    "$CONTRASTS_FILE" \
    "$GMT_FILE"

# --- STEP 5: FINAL MULTIQC AGGREGATION ---
echo "RUNNING STEP 5: FINAL MULTIQC AGGREGATION"

cat << 'CONFIG' > mqc_config.yaml
# mqc_config.yaml
disable_version_detection: true
sections:
  software_versions:
    hide: true

custom_content:
  order:
    - pathway_analysis
CONFIG

MULTIQC_CONTAINER="docker://multiqc/multiqc:v1.33"

singularity exec --bind $PWD:/data --pwd /data $MULTIQC_CONTAINER multiqc \
    /data/ANALYSIS/results_environment \
    /data/ANALYSIS/results_human_final \
    /data/ANALYSIS/xengsort_out \
    --force \
    --config /data/mqc_config.yaml \
    --title "U251 Environment Comparison" \
    --filename "U251_Environment_Report.html" \
    --outdir "/data/ANALYSIS/results_environment"

echo "DONE."
