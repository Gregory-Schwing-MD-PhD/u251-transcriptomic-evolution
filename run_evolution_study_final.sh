#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=go2432@wayne.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o xengsort_prod_%j.out
#SBATCH -e xengsort_prod_%j.err

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

# --- STEP 1: XENGSORT (Deconvolution) ---
# UPDATED: Using params file for clean execution
echo "RUNNING STEP 1: XENGSORT"
nextflow run main.nf \
    -profile singularity \
    -params-file xengsort_params.yaml \
    -w "${WORK_DIR}" \
    -resume \
    -ansi-log false

unset NXF_PARAMS

# --- STEP 2: RNA-SEQ (Quantification & QC) ---
echo "RUNNING STEP 2: RNA-SEQ"
nextflow run nf-core/rnaseq \
    -r 3.22.2 \
    -profile singularity \
    -c rnaseq_custom.config \
    -params-file rnaseq_params.yaml \
    -w "${WORK_DIR}" \
    -resume \
    -ansi-log false

unset NXF_PARAMS

# --- STEP 3: DIFFERENTIAL ABUNDANCE ---
echo "RUNNING STEP 3: DIFFERENTIAL ABUNDANCE"
nextflow run nf-core/differentialabundance \
    -r 1.5.0 \
    -profile singularity \
    --input "$(pwd)/ANALYSIS/metadata_therapy.csv" \
    --contrasts "$(pwd)/ANALYSIS/contrasts_therapy.csv" \
    --matrix "$(pwd)/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv" \
    --transcript_length_matrix "$(pwd)/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_lengths.tsv" \
    --gtf "$(pwd)/ANALYSIS/refs/human/GRCh38.primary_assembly.annotation.gtf.gz" \
    --exploratory_main_variable Classification \
    --outdir "ANALYSIS/results_therapy" \
    -params-file therapy_params.yaml \
    -w "${WORK_DIR}" \
    -resume \
    -ansi-log false

unset NXF_PARAMS

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
