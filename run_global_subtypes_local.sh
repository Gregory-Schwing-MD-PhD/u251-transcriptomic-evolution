#!/bin/bash
set -euo pipefail

# ==============================================================================
# U251 EVOLUTIONARY ANALYSIS - ULTIMATE EDITION v16.4
# ==============================================================================

echo "========================================"
echo "U251 Transcriptomic Evolution Analysis"
echo "Ultimate Edition v16.4"
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
echo "Running Ultimate Analysis v16.4"
echo "Features:"
echo "  • 2×4×N Statistical Matrix"
echo "  • Comprehensive Pairwise Testing"
echo "  • Enhanced LLM Summary (embedded)"
echo "  • Professional Significance Heatmap"
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
title: "U251 Global Evolution: Ultimate Analysis v16.4"
subtitle: "Comprehensive Transcriptomic Evolution with 2×4×N Statistical Matrix"
intro_text: |
  This report presents the ultimate comprehensive analysis of U251 glioblastoma 
  transcriptomic evolution across Culture → Primary → Recurrent stages.

  **Key Features:**
  - 2×4×N Statistical Framework (GSVA/Z-Score × 4 Tests × All Comparisons)
  - Professional significance matrix with visual heatmap
  - Complete correlation analysis (all signature pairs)
  - Enhanced LLM interpretation summary (embedded in report)
  - PCA gene drivers and plasticity statistics
  - Jonckheere-Terpstra trajectory testing with consensus

  **Statistical Rigor:** Multiple testing correction (FDR), dual scoring methods,
  parametric and non-parametric tests, weighted and unweighted approaches.

report_header_info:
  - Pipeline Version: "v16.4 Ultimate Edition"
  - Statistical Framework: "2×4×N Matrix"
  - Scoring Methods: "GSVA + Z-Score"
  - Trajectory Testing: "Jonckheere-Terpstra + Polynomial"

custom_content:
  order:
    - analysis_summary
    - significance_matrix
    - trajectory_tests
    - plasticity_matrix
    - global_structure
    - method_agreement

extra_fn_clean_exts:
  - "_mqc"

table_columns_visible:
  significance_matrix:
    all_columns: true
  plasticity_matrix:
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
    --filename "U251_Ultimate_Analysis_v16.4.html" \
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
echo "ULTIMATE ANALYSIS v16.4 COMPLETE"
echo "========================================"
echo ""
echo "Primary Output:"
echo "  📊 Report: ${OUT_DIR}/U251_Ultimate_Analysis_v16.4.html"
echo "  📄 LLM Summary: ${OUTPUT_PREFIX}_llm_summary_ULTIMATE.txt"
echo ""
echo "Key Data Files:"
echo "  • Significance Matrix (Full): ${OUTPUT_PREFIX}_Significance_Matrix_Full.csv"
echo "  • Trajectory Tests (All): ${OUTPUT_PREFIX}_Trajectory_Comprehensive_ALL_Tests.csv"
echo "  • Plasticity Matrix (2×4): ${OUTPUT_PREFIX}_Plasticity_Comprehensive_Matrix.csv"
echo "  • Correlation Matrix (Complete): ${OUTPUT_PREFIX}_Signature_Correlation_Matrix_FULL.csv"
echo ""
echo "Visualizations:"
ls -1 ${OUTPUT_PREFIX}*_mqc.png 2>/dev/null | while read file; do
    echo "  • $(basename $file)"
done
echo ""
echo "Session Info: ${OUT_DIR}/sessionInfo.txt"
echo ""
echo "========================================"
echo "Next Steps:"
echo "  1. Open HTML report in browser"
echo "  2. Review embedded LLM summary section"
echo "  3. Copy LLM summary to ChatGPT/Claude for interpretation"
echo "  4. Examine significance matrix heatmap"
echo "  5. Check trajectory consensus trends"
echo "========================================"
date
echo "Analysis completed successfully!"
echo "========================================"
