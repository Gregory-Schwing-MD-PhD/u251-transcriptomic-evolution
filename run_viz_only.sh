#!/bin/bash
#SBATCH -q secondary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:45:00
#SBATCH --job-name=u251_viz_ks
#SBATCH -o viz_ks_%j.out
#SBATCH -e viz_ks_%j.err
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
# 4. DEFINE PATHS
# ==========================================
CONTAINER="docker://go2432/bioconductor:latest"
IMG_PATH="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"
MULTIQC_CONTAINER="docker://multiqc/multiqc:v1.33"

# Inputs from your previous successful run
RESULTS_DIR="ANALYSIS/results_evolution"
OUTPUT_PREFIX="ANALYSIS/results_evolution/plots/Evolution"
COUNTS_FILE="ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv"
METADATA_FILE="ANALYSIS/metadata.csv"
CONTRASTS_FILE="ANALYSIS/contrasts.csv"
GMT_FILE="ANALYSIS/refs/pathways/combined_human.gmt"

# Pull container if missing (safety check)
if [[ ! -f "$IMG_PATH" ]]; then
    singularity pull "$IMG_PATH" "$CONTAINER"
fi

# ==========================================
# 5. RUN THE EVOLUTIONARY KITCHEN SINK (R)
# ==========================================
echo "RUNNING R VISUALIZATION SUITE..."

singularity exec --bind $PWD:/data --pwd /data "$IMG_PATH" Rscript plot_evolution_kitchen_sink.R \
    "$RESULTS_DIR" \
    "$OUTPUT_PREFIX" \
    "$COUNTS_FILE" \
    "$METADATA_FILE" \
    "$CONTRASTS_FILE" \
    "$GMT_FILE"

# ==========================================
# 6. RUN MULTIQC
# ==========================================
echo "RUNNING FINAL MULTIQC AGGREGATION..."

cat << 'CONFIG' > mqc_config.yaml
disable_version_detection: true
sections:
  software_versions:
    hide: true
custom_content:
  order:
    - pathway_analysis
CONFIG

singularity exec --bind $PWD:/data --pwd /data $MULTIQC_CONTAINER multiqc \
    /data/ANALYSIS/results_evolution \
    /data/ANALYSIS/results_human_final \
    /data/ANALYSIS/xengsort_out \
    --force \
    --config /data/mqc_config.yaml \
    --title "U251 Transcriptomic Evolution" \
    --filename "U251_Evolution_Report.html" \
    --outdir "/data/ANALYSIS/results_evolution"

echo "DONE."
