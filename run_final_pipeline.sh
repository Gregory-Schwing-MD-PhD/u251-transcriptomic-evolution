#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=2:00:00
#SBATCH --job-name=u251_final
#SBATCH --output=final_analysis_%j.out
#SBATCH --error=final_analysis_%j.err

set -x

# 1. SETUP (Copied exact environment config)
export CONDA_PREFIX="${HOME}/mambaforge/envs/nextflow"
export PATH="${CONDA_PREFIX}/bin:$PATH"
unset JAVA_HOME
nextflow -v

# 2. CACHE
export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
mkdir -p $XDG_RUNTIME_DIR $NXF_SINGULARITY_CACHEDIR

# 3. SAFETY
export NXF_SINGULARITY_HOME_MOUNT=true
unset LD_LIBRARY_PATH
unset PYTHONPATH
unset R_LIBS
unset R_LIBS_USER
unset R_LIBS_SITE

# ==============================================================================
# 4. PATH CONFIGURATION
# ==============================================================================
GMT_DIR="ANALYSIS/refs/pathways/human"
STRING_DIR="ANALYSIS/refs/pathways/human/human_string"
RESULTS_DIR="ANALYSIS/results_evolution"
VST_FILE="${RESULTS_DIR}/tables/processed_abundance/all.vst.tsv"

# The Clean Output Directory for MultiQC
OUT_DIR="ANALYSIS/results_visualization/Final_Report"
mkdir -p "$OUT_DIR"

# Singularity Container
IMG_PATH="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"
MULTIQC_CONTAINER="docker://multiqc/multiqc:v1.33"

if [[ ! -f "$IMG_PATH" ]]; then
    singularity pull "$IMG_PATH" docker://go2432/bioconductor:latest
fi

# ==============================================================================
# 5. EXECUTE R ANALYSIS (Generates _mqc files)
# ==============================================================================
echo "RUNNING STEP 1: PATHWAYS, DRUGS, & STRING PPI"

singularity exec --bind $PWD:/data --pwd /data "$IMG_PATH" \
    Rscript run_pathways_drugs_v4.R \
    "$VST_FILE" \
    "$RESULTS_DIR" \
    "$GMT_DIR" \
    "$STRING_DIR" \
    "$OUT_DIR/Analysis"

# ==============================================================================
# 6. GENERATE MULTIQC REPORT
# ==============================================================================
echo "RUNNING STEP 2: MULTIQC AGGREGATION"

# Create Custom Config to prioritize our custom content
cat << 'CONFIG' > mqc_config.yaml
disable_version_detection: true
sections:
  software_versions:
    hide: true
custom_content:
  order:
    - analysis_narrative
    - pathway_analysis
sp:
  analysis_narrative:
    fn: "*_mqc.html"
CONFIG

# Run MultiQC
# It will scan $OUT_DIR, find Analysis_Narrative_mqc.html, and include it.
singularity exec --bind $PWD:/data --pwd /data $MULTIQC_CONTAINER multiqc \
    /data/$OUT_DIR \
    --force \
    --config /data/mqc_config.yaml \
    --title "U251 Transcriptomic Evolution & Drug Discovery" \
    --filename "U251_Final_Report.html" \
    --outdir "/data/$OUT_DIR"

echo "--------------------------------------------------------"
echo "ANALYSIS COMPLETE"
echo "Report: $OUT_DIR/U251_Final_Report.html"
echo "--------------------------------------------------------"
