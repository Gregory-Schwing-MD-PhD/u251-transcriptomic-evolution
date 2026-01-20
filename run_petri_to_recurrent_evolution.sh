#!/bin/bash
#SBATCH -q secondary
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

# 4. EXECUTION
sed -i 's/\r//g' ANALYSIS/metadata.csv

echo "RUNNING STEP 3: DIFFERENTIAL ABUNDANCE (EVOLUTION)"
# The Pipeline does all the math here!
nextflow run nf-core/differentialabundance \
    -r 1.5.0 \
    -profile singularity \
    -params-file evolution_params.yaml \
    -w "${WORK_DIR}" \
    -ansi-log false \
    -resume

# 5. VISUALIZATION (PURE VIZ)
echo "RUNNING STEP 4: KITCHEN SINK (VISUALIZATION ONLY)"

CONTAINER="docker://go2432/bioconductor:latest"
IMG_PATH="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"

RESULTS_DIR="ANALYSIS/results_evolution"
OUTPUT_PREFIX="ANALYSIS/results_evolution/plots/Evolution"
METADATA_FILE="ANALYSIS/metadata.csv"
CONTRASTS_FILE="ANALYSIS/contrasts.csv"
GMT_FILE="ANALYSIS/refs/pathways/combined_human.gmt"

if [[ ! -f "$IMG_PATH" ]]; then
    singularity pull "$IMG_PATH" "$CONTAINER"
fi

# R Arguments: <results_dir> <out_prefix> <meta> <contrasts> <gmt>
singularity exec --bind $PWD:/data --pwd /data "$IMG_PATH" Rscript plot_evolution_kitchen_sink.R \
    "$RESULTS_DIR" \
    "$OUTPUT_PREFIX" \
    "$METADATA_FILE" \
    "$CONTRASTS_FILE" \
    "$GMT_FILE"

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
