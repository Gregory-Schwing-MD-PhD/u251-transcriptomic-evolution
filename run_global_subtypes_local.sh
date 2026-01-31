#!/bin/bash
set -euo pipefail

# ==============================================================================
# U251 EVOLUTIONARY ANALYSIS WITH MULTI-FDR SENSITIVITY
# ==============================================================================

echo "========================================"
echo "U251 Transcriptomic Evolution Analysis"
echo "with Multi-FDR Sensitivity Testing"
echo "========================================"
date

# ------------------------------------------------------------------------------
# ENVIRONMENT SETUP
# ------------------------------------------------------------------------------
export CONDA_PREFIX="${HOME}/mambaforge/envs/nextflow"
export PATH="${CONDA_PREFIX}/bin:$PATH"
export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
mkdir -p $XDG_RUNTIME_DIR $NXF_SINGULARITY_CACHEDIR
export NXF_SINGULARITY_HOME_MOUNT=true
unset LD_LIBRARY_PATH
unset R_LIBS

# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------
CONTAINER="docker://go2432/bioconductor:latest"
IMG_PATH="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"

RESULTS_DIR="ANALYSIS/results_evolution"
METADATA_FILE="ANALYSIS/metadata.csv"
OUT_DIR="ANALYSIS/results_global_subtypes"
OUTPUT_PREFIX="${OUT_DIR}/Global_Subtypes"

# ------------------------------------------------------------------------------
# INPUT VALIDATION
# ------------------------------------------------------------------------------
echo ""
echo "Validating inputs..."

VST_FILE="${RESULTS_DIR}/tables/processed_abundance/all.vst.tsv"

for FILE in "$METADATA_FILE" "$VST_FILE"; do
    if [[ ! -f "$FILE" ]]; then
        echo "ERROR: Required file not found: $FILE"
        exit 1
    fi
    echo "  ✓ Found: $FILE"
done

if [[ ! -f "run_global_subtypes.R" ]]; then
    echo "ERROR: run_global_subtypes.R not found in current directory"
    exit 1
fi
echo "  ✓ Found: run_global_subtypes.R"

mkdir -p "$OUT_DIR"
echo "  ✓ Output directory: $OUT_DIR"

# ------------------------------------------------------------------------------
# CONTAINER MANAGEMENT
# ------------------------------------------------------------------------------
echo ""
echo "Managing container..."

if [[ ! -f "$IMG_PATH" ]]; then
    echo "  Pulling container (this may take a few minutes)..."
    singularity pull "$IMG_PATH" "$CONTAINER" || {
        echo "ERROR: Failed to pull container"
        exit 1
    }
    echo "  ✓ Container pulled successfully"
else
    echo "  ✓ Using cached container: $IMG_PATH"
fi

# ------------------------------------------------------------------------------
# RUN R ANALYSIS
# ------------------------------------------------------------------------------
echo ""
echo "========================================"
echo "Running Multi-FDR Sensitivity Analysis"
echo "Testing at α = 0.05, 0.01, 0.005, 0.001"
echo "========================================"

singularity exec --bind $PWD:/data --pwd /data "$IMG_PATH" \
    Rscript run_global_subtypes.R \
    --input "$VST_FILE" \
    --out "$OUTPUT_PREFIX" \
    --meta "$METADATA_FILE" || {
        echo "ERROR: R analysis failed"
        exit 1
    }

echo ""
echo "✓ R analysis completed successfully"

# ------------------------------------------------------------------------------
# GENERATE MULTIQC REPORT
# ------------------------------------------------------------------------------
echo ""
echo "========================================"
echo "Generating MultiQC Report..."
echo "========================================"

cat > mqc_config_subtypes.yaml << 'CONFIG_EOF'
title: "U251 Global Evolution: Multi-FDR Sensitivity Report"
subtitle: "Transcriptomic Evolution with Multiple Statistical Stringency Levels"
intro_text: |
  This report presents a comprehensive multi-FDR sensitivity analysis of U251 
  glioblastoma transcriptomic evolution. Due to small sample sizes (n=1-3 per group), 
  results are shown at four FDR thresholds (α=0.05, 0.01, 0.005, 0.001) to balance 
  discovery and false positive control.
  
  **Recommended Interpretation:** Focus on signatures significant at α≤0.01 and 
  validate findings with orthogonal experimental methods.

report_header_info:
  - Pipeline Version: "v12.0 (Multi-FDR)"
  - FDR Thresholds Tested: "0.05, 0.01, 0.005, 0.001"

custom_content:
  order:
    - fdr_sensitivity
    - subtype_stats_fdr_0001
    - subtype_stats_fdr_0005  
    - subtype_stats_fdr_001
    - subtype_stats_fdr_005
    - analysis_summary

extra_fn_clean_exts:
  - "_mqc"

table_columns_visible:
  fdr_sensitivity:
    all_columns: true
  subtype_stats_fdr_005:
    all_columns: true

CONFIG_EOF

MULTIQC_CONTAINER="docker://multiqc/multiqc:v1.33"
MULTIQC_IMG="${NXF_SINGULARITY_CACHEDIR}/multiqc.sif"

if [[ ! -f "$MULTIQC_IMG" ]]; then
    echo "  Pulling MultiQC container..."
    singularity pull "$MULTIQC_IMG" "$MULTIQC_CONTAINER" || {
        echo "WARNING: Failed to pull MultiQC container"
    }
fi

singularity exec --bind $PWD:/data --pwd /data "$MULTIQC_IMG" multiqc \
    "$OUT_DIR" \
    --force \
    --config mqc_config_subtypes.yaml \
    --filename "U251_Multi_FDR_Sensitivity_Report.html" \
    --outdir "$OUT_DIR" || {
        echo "ERROR: MultiQC failed"
        exit 1
    }

echo ""
echo "✓ MultiQC report generated successfully"

# ------------------------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------------------------
echo ""
echo "========================================"
echo "MULTI-FDR ANALYSIS COMPLETE"
echo "========================================"
echo ""
echo "Primary Output:"
echo "  📊 Report: ${OUT_DIR}/U251_Multi_FDR_Sensitivity_Report.html"
echo ""
echo "Statistics Tables (by FDR threshold):"
ls -1 ${OUT_DIR}/fdr_*_stats_mqc.tsv 2>/dev/null | while read file; do
    echo "  - $(basename $file)"
done
echo ""
echo "Data Exports:"
ls -1 ${OUTPUT_PREFIX}*.csv 2>/dev/null | while read file; do
    echo "  - $(basename $file)"
done
echo ""
echo "Figures:"
ls -1 ${OUTPUT_PREFIX}*.png 2>/dev/null | while read file; do
    echo "  - $(basename $file)"
done
echo ""
echo "Session Info: ${OUT_DIR}/sessionInfo.txt"
echo ""
echo "========================================"
date
echo "Analysis completed successfully!"
echo "========================================"
