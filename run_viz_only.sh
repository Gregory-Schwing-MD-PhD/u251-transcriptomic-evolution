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
# 4. VISUALIZATION (THE KITCHEN SINK)
# ==========================================
echo "RUNNING STEP 4: GENERATING PUBLICATION PLOTS"

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
COUNTS_FILE="ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv"

# FIX: Use $PWD:/data to provide absolute paths and set the working directory
# This resolves the "destination must be an absolute path" FATAL error
singularity exec --bind $PWD:/data --pwd /data "$IMG_PATH" Rscript plot_kitchen_sink.R \
    "$DESEQ_FILE" \
    "$VST_FILE" \
    "$GMT_FILE" \
    "$OUTPUT_PREFIX" \
    "$COUNTS_FILE"

# ==========================================
# 5. FINAL MULTIQC AGGREGATION
# ==========================================
echo "RUNNING STEP 5: FINAL MULTIQC AGGREGATION"

# FIX: MultiQC v1.33 KeyError bypass
# We create a local config to explicitly disable the software_versions module 
# before the Python Traceback occurs.
cat << 'CONFIG' > mqc_config.yaml
software_versions:
    show: false
CONFIG

MULTIQC_CONTAINER="docker://multiqc/multiqc:v1.33"

# Run MultiQC using absolute paths
# We point to the specific directories to ensure it finds the data in the mapped /data volume
singularity exec --bind $PWD:/data --pwd /data $MULTIQC_CONTAINER multiqc \
    /data/ANALYSIS/results_therapy \
    /data/ANALYSIS/results_human_final \
    /data/ANALYSIS/xengsort_out \
    --force \
    --config /data/mqc_config.yaml \
    --title "U251 Transcriptomic Evolution" \
    --filename "U251_Final_Report.html" \
    --outdir "/data/ANALYSIS/results_therapy" \
    --exclude software_versions

echo "DONE."
