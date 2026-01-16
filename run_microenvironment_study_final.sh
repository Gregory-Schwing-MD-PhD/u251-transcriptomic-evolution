#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=go2432@wayne.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o host_analysis_%j.out
#SBATCH -e host_analysis_%j.err

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

# Independent work directory for the host analysis
WORK_DIR="$(pwd)/work_host"

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

# --- STEP 1: XENGSORT (SKIPPED) ---
# echo "SKIPPING STEP 1: XENGSORT"
# nextflow run main.nf \
#     -profile singularity \
#     -params-file xengsort_params.yaml \
#     -w "${WORK_DIR}" \
#     -resume \
#     -ansi-log false

# unset NXF_PARAMS

# --- STEP 2: HOST RNA-SEQ (ACTIVE) ---
# Processes the Rat-specific FASTQs
echo "RUNNING STEP 2: HOST RNA-SEQ"
nextflow run nf-core/rnaseq \
    -r 3.22.2 \
    -profile singularity \
    -c rnaseq_custom.config \
    -params-file host_params.yaml \
    -w "${WORK_DIR}" \
    -resume \
    -ansi-log false

unset NXF_PARAMS

# --- STEP 3: DIFFERENTIAL ABUNDANCE ---
echo "RUNNING STEP 3: DIFFERENTIAL ABUNDANCE"
nextflow run nf-core/differentialabundance \
    -r 1.5.0 \
    -profile singularity \
    --input "$(pwd)/ANALYSIS/metadata_host.csv" \
    --contrasts "$(pwd)/ANALYSIS/contrasts_host.csv" \
    --matrix "$(pwd)/ANALYSIS/results_host_rat/star_salmon/salmon.merged.gene_counts.tsv" \
    --transcript_length_matrix "$(pwd)/ANALYSIS/results_host_rat/star_salmon/salmon.merged.gene_lengths.tsv" \
    --gtf "$(pwd)/ANALYSIS/refs/rat/Rattus_norvegicus.mRatBN7.2.110.gtf.gz" \
    --exploratory_main_variable condition \
    --outdir "ANALYSIS/results_host_rat/differential" \
    -params-file therapy_params.yaml \
    -w "${WORK_DIR}" \
    -resume \
    -ansi-log false

unset NXF_PARAMS

# --- STEP 4: SIGNAL SUBTRACTION (HOST SPECIFIC) ---
echo "RUNNING STEP 4: SIGNAL SUBTRACTION"

CONTAINER="docker://go2432/bioconductor:latest"
IMG_PATH="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"

# Pull if not exists
if [[ ! -f "$IMG_PATH" ]]; then
    singularity pull "$IMG_PATH" "$CONTAINER"
fi

# We use the same container setup as the visualization step
singularity exec --bind $PWD:/data --pwd /data "$IMG_PATH" Rscript plot_host_cleaner.R \
    "ANALYSIS/results_host_rat/star_salmon/salmon.merged.gene_counts.tsv" \
    "ANALYSIS/metadata_host.csv" \
    "ANALYSIS/results_host_rat/LITT_Microenvironment"

# --- STEP 5: VISUALIZATION (THE KITCHEN SINK) ---
echo "RUNNING STEP 5: GENERATING PUBLICATION PLOTS"

# Inputs pointing to RAT results
DESEQ_FILE="ANALYSIS/results_host_rat/differential/tables/differential/LITT_Response.deseq2.results.tsv"
VST_FILE="ANALYSIS/results_host_rat/differential/tables/processed_abundance/all.vst.tsv"
GMT_FILE="ANALYSIS/refs/pathways/combined_rat.gmt"
OUTPUT_PREFIX="ANALYSIS/results_host_rat/plots/Rat_Microenvironment"
COUNTS_FILE="ANALYSIS/results_host_rat/star_salmon/salmon.merged.gene_counts.tsv"

# NEW INPUT: The clean gene list generated in Step 4
CLEAN_GENES_FILE="ANALYSIS/results_host_rat/LITT_Microenvironment_Clean_LITT_Genes.csv"

# Verification
if [[ ! -f "$DESEQ_FILE" ]]; then
    echo "ERROR: DESeq2 results file not found at: $DESEQ_FILE"
    exit 1
fi

# UPDATED: Calling 'plot_kitchen_sink_microenvironment.R' with the 6th argument
singularity exec --bind $PWD:/data --pwd /data "$IMG_PATH" Rscript plot_kitchen_sink_microenvironment.R \
    "$DESEQ_FILE" \
    "$VST_FILE" \
    "$GMT_FILE" \
    "$OUTPUT_PREFIX" \
    "$COUNTS_FILE" \
    "$CLEAN_GENES_FILE"

# --- STEP 6: GENERATE FINAL MULTIQC REPORT ---
echo "RUNNING STEP 6: FINAL MULTIQC AGGREGATION"

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

# Run MultiQC using absolute paths via /data bind
# Removed --exclude software_versions flag as we handle it in config now
singularity exec --bind $PWD:/data --pwd /data $MULTIQC_CONTAINER multiqc \
    /data/ANALYSIS/results_host_rat \
    /data/ANALYSIS/xengsort_out \
    --force \
    --config /data/mqc_config.yaml \
    --title "U251 Microenvironment Response (Rat)" \
    --filename "U251_Host_Report.html" \
    --outdir "/data/ANALYSIS/results_host_rat"

echo "DONE. Final report: ANALYSIS/results_host_rat/U251_Host_Report.html"
