#!/bin/bash
#SBATCH -q secondary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH -o viz_only_%j.out
#SBATCH -e viz_only_%j.err

set -x

# ==========================================
# 1. ENVIRONMENT SETUP
# ==========================================
export CONDA_PREFIX="${HOME}/mambaforge/envs/nextflow"
export PATH="${CONDA_PREFIX}/bin:$PATH"
unset JAVA_HOME

export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
mkdir -p $XDG_RUNTIME_DIR $NXF_SINGULARITY_CACHEDIR

export NXF_SINGULARITY_HOME_MOUNT=true
unset LD_LIBRARY_PATH PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE

# ==========================================
# 2. VISUALIZATION (THE KITCHEN SINK)
# ==========================================
echo "RUNNING STEP 2: GENERATING PUBLICATION PLOTS"

CONTAINER="docker://go2432/bioconductor:latest"
IMG_PATH="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"

# Pull if not exists
if [[ ! -f "$IMG_PATH" ]]; then
    singularity pull "$IMG_PATH" "$CONTAINER"
fi

# Inputs
DESEQ_FILE="ANALYSIS/results_therapy/tables/differential/therapy_impact.deseq2.results.tsv"
VST_FILE="ANALYSIS/results_therapy/tables/processed_abundance/all.vst.tsv"
GMT_FILE="ANALYSIS/refs/pathways/combined_human.gmt"
OUTPUT_PREFIX="ANALYSIS/results_therapy/plots/U251_Publication"

# UPDATED: Removed counts_file argument (v5 script uses internal org.Hs.eg.db mapping)
singularity exec --bind $PWD:/data --pwd /data "$IMG_PATH" Rscript plot_kitchen_sink.R \
    "$DESEQ_FILE" \
    "$VST_FILE" \
    "$GMT_FILE" \
    "$OUTPUT_PREFIX"

# ==========================================
# 3. FINAL MULTIQC AGGREGATION
# ==========================================
echo "RUNNING STEP 3: FINAL MULTIQC AGGREGATION"

# Config to silence software_versions internally
cat << 'CONFIG' > mqc_config.yaml
# mqc_config.yaml
disable_version_detection: true
sections:
  software_versions:
    hide: true

run_modules:
  - custom_content

custom_content:
  order:
    - pathway_analysis

CONFIG

MULTIQC_CONTAINER="docker://multiqc/multiqc:v1.33"

# Run MultiQC
singularity exec --bind $PWD:/data --pwd /data $MULTIQC_CONTAINER multiqc \
    /data/ANALYSIS/results_therapy \
    /data/ANALYSIS/results_human_final \
    /data/ANALYSIS/xengsort_out \
    --force \
    --config /data/mqc_config.yaml \
    --title "U251 Transcriptomic Evolution" \
    --filename "U251_Final_Report.html" \
    --outdir "/data/ANALYSIS/results_therapy"

echo "DONE."
