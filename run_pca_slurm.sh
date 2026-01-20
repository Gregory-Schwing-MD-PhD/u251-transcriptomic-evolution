#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --job-name=u251_pca_pub
#SBATCH -o pca_pub_%j.out
#SBATCH -e pca_pub_%j.err
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

# ==========================================
# 2. CACHE & STORAGE SETUP
# ==========================================
export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
mkdir -p $XDG_RUNTIME_DIR $NXF_SINGULARITY_CACHEDIR

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
# 4. DEFINE PATHS
# ==========================================
CONTAINER="docker://go2432/bioconductor:latest"
IMG_PATH="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"
MULTIQC_CONTAINER="docker://multiqc/multiqc:v1.33"

# Inputs
COUNTS_FILE="ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv"
METADATA_FILE="ANALYSIS/metadata.csv"
GMT_FILE="ANALYSIS/refs/pathways/combined_human.gmt"

# Output
RESULTS_DIR="ANALYSIS/results_pca_pub"
OUTPUT_PREFIX="${RESULTS_DIR}/plots/U251"

if [[ ! -f "$IMG_PATH" ]]; then
    singularity pull "$IMG_PATH" "$CONTAINER"
fi

# ==========================================
# 5. RUN CUSTOM PCA.R
# ==========================================
echo "RUNNING PCA PUBLICATION SUITE..."

singularity exec --bind $PWD:/data --pwd /data "$IMG_PATH" Rscript PCA.R \
    "$COUNTS_FILE" \
    "$METADATA_FILE" \
    "$GMT_FILE" \
    "$OUTPUT_PREFIX"

# ==========================================
# 6. RUN MULTIQC
# ==========================================
echo "RUNNING MULTIQC AGGREGATION..."

singularity exec --bind $PWD:/data --pwd /data $MULTIQC_CONTAINER multiqc \
    "/data/${RESULTS_DIR}" \
    --force \
    --config "/data/multiqc_config_pca.yaml" \
    --title "U251 PCA & Trajectory" \
    --filename "U251_PCA_Report.html" \
    --outdir "/data/${RESULTS_DIR}"

echo "DONE."
