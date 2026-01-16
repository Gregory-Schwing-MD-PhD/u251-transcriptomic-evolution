#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --job-name=rat_evolution
#SBATCH -o host_analysis_%j.out
#SBATCH -e host_analysis_%j.err

set -x

# SETUP
export CONDA_PREFIX="${HOME}/mambaforge/envs/nextflow"
export PATH="${CONDA_PREFIX}/bin:$PATH"
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
WORK_DIR="$(pwd)/work_host"

# --- STEP 4: SIGNAL SUBTRACTION (Ensure this runs first to get clean genes) ---
echo "RUNNING STEP 4: SIGNAL SUBTRACTION"

CONTAINER="docker://go2432/bioconductor:latest"
IMG_PATH="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"

if [[ ! -f "$IMG_PATH" ]]; then
    singularity pull "$IMG_PATH" "$CONTAINER"
fi

# Note: This script (plot_host_cleaner.R) generates the Clean Genes list
singularity exec --bind $PWD:/data --pwd /data "$IMG_PATH" Rscript plot_host_cleaner.R \
    "ANALYSIS/results_host_rat/star_salmon/salmon.merged.gene_counts.tsv" \
    "ANALYSIS/metadata_host.csv" \
    "ANALYSIS/results_host_rat/LITT_Microenvironment"

# --- STEP 5: VISUALIZATION (RAT EVOLUTION) ---
echo "RUNNING STEP 5: RAT EVOLUTION VISUALIZATION"

# Inputs
RESULTS_DIR="ANALYSIS/results_host_rat"
OUTPUT_PREFIX="ANALYSIS/results_host_rat/plots/Rat_Evolution"
COUNTS_FILE="ANALYSIS/results_host_rat/star_salmon/salmon.merged.gene_counts.tsv"
METADATA_FILE="ANALYSIS/metadata_host.csv"
CONTRASTS_FILE="ANALYSIS/contrasts_host.csv"
GMT_FILE="ANALYSIS/refs/pathways/combined_rat.gmt" # Using RAT GMT
CLEAN_GENES="ANALYSIS/results_host_rat/LITT_Microenvironment_Clean_LITT_Genes.csv"

# Run the NEW Rat Evolution Script
singularity exec --bind $PWD:/data --pwd /data "$IMG_PATH" Rscript plot_host_evolution.R \
    "$RESULTS_DIR" \
    "$OUTPUT_PREFIX" \
    "$COUNTS_FILE" \
    "$METADATA_FILE" \
    "$CONTRASTS_FILE" \
    "$GMT_FILE" \
    "$CLEAN_GENES"

# --- STEP 6: FINAL MULTIQC AGGREGATION ---
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

# Run MultiQC
singularity exec --bind $PWD:/data --pwd /data $MULTIQC_CONTAINER multiqc \
    /data/ANALYSIS/results_host_rat \
    /data/ANALYSIS/xengsort_out \
    --force \
    --config /data/mqc_config.yaml \
    --title "U251 Microenvironment Response (Rat)" \
    --filename "U251_Host_Report.html" \
    --outdir "/data/ANALYSIS/results_host_rat"

echo "DONE."
