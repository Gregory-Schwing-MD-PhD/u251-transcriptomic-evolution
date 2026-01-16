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

set -x

# SETUP
export CONDA_PREFIX="${HOME}/mambaforge/envs/nextflow"
export PATH="${CONDA_PREFIX}/bin:$PATH"
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
WORK_DIR="$(pwd)/work"

# --- SKIP STEPS 1 & 2 (Using existing results) ---
echo "LOG: Skipping Xengsort and RNA-Seq Alignment. Using existing matrices."

# --- STEP 3: DIFFERENTIAL ABUNDANCE (EVOLUTION RUN) ---
echo "RUNNING STEP 3: DIFFERENTIAL ABUNDANCE (EVOLUTION)"

nextflow run nf-core/differentialabundance \
    -r 1.5.0 \
    -profile singularity \
    --input "$(pwd)/ANALYSIS/metadata.csv" \
    --contrasts "$(pwd)/ANALYSIS/contrasts.csv" \
    --matrix "$(pwd)/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv" \
    --transcript_length_matrix "$(pwd)/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_lengths.tsv" \
    --gtf "$(pwd)/ANALYSIS/refs/human/GRCh38.primary_assembly.annotation.gtf.gz" \
    --exploratory_main_variable Classification \
    --outdir "ANALYSIS/results_evolution" \
    -w "${WORK_DIR}" \
    -resume \
    -ansi-log false

# --- STEP 4: CUSTOM VISUALIZATION (MEGA SCRIPT) ---
echo "RUNNING STEP 4: KITCHEN SINK + EVOLUTION"

CONTAINER="docker://go2432/bioconductor:latest"
IMG_PATH="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"

# Inputs for R script
RESULTS_DIR="ANALYSIS/results_evolution"
OUTPUT_PREFIX="ANALYSIS/results_evolution/plots/Evolution"
COUNTS_FILE="ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv"
METADATA_FILE="ANALYSIS/metadata.csv"
CONTRASTS_FILE="ANALYSIS/contrasts.csv"
GMT_FILE="ANALYSIS/refs/pathways/combined_human.gmt"

if [[ ! -f "$IMG_PATH" ]]; then
    singularity pull "$IMG_PATH" "$CONTAINER"
fi

# Pass ALL 6 arguments (including GMT)
singularity exec --bind $PWD:/data --pwd /data "$IMG_PATH" Rscript plot_evolution_kitchen_sink.R \
    "$RESULTS_DIR" \
    "$OUTPUT_PREFIX" \
    "$COUNTS_FILE" \
    "$METADATA_FILE" \
    "$CONTRASTS_FILE" \
    "$GMT_FILE"

# --- STEP 5: FINAL MULTIQC AGGREGATION ---
echo "RUNNING STEP 5: FINAL MULTIQC AGGREGATION"

# UPDATED: Corrected Config Structure
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

# UPDATED: 
# 1. Scans 'results_evolution' (where Step 4 put the plots)
# 2. outputs report to 'results_evolution'
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
