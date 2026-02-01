#!/bin/bash

set -euo pipefail

# ==============================================================================
# PUBLICATION FIGURE GENERATION SCRIPT
# ==============================================================================
# Generates 9-panel publication figure (A-I) with comprehensive visualizations
# ==============================================================================

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║      PUBLICATION FIGURE GENERATOR - 9 PANEL COMPOSITE          ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

# ==============================================================================
# SETUP
# ==============================================================================
export CONDA_PREFIX="${HOME}/mambaforge/envs/nextflow"
export PATH="${CONDA_PREFIX}/bin:$PATH"
unset JAVA_HOME

# Cache directories
export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
mkdir -p "$XDG_RUNTIME_DIR" "$NXF_SINGULARITY_CACHEDIR"

# Clean environment
export NXF_SINGULARITY_HOME_MOUNT=true
unset LD_LIBRARY_PATH PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE

# ==============================================================================
# PATH CONFIGURATION
# ==============================================================================
# CORRECTED PATH TO COMBINED GMT
GMT_DIR="ANALYSIS/refs/pathways/human"
COMBINED_GMT="${GMT_DIR}/combined_human.gmt"

STRING_DIR="ANALYSIS/refs/pathways/human/human_string"
RESULTS_DIR="ANALYSIS/results_evolution"
VST_FILE="${RESULTS_DIR}/tables/processed_abundance/all.vst.tsv"
METADATA_FILE="ANALYSIS/metadata.csv"  # ADDED METADATA FILE
OUT_DIR="publication_figure"
TARGET_CONTRAST="therapy_impact"

# Verify paths
echo "Verifying input files..."
for path in "$VST_FILE" "$COMBINED_GMT" "$STRING_DIR" "$METADATA_FILE"; do
    if [ ! -e "$path" ]; then
        echo "ERROR: Required path not found: $path"
        exit 1
    fi
done

if [ ! -f "Experiment_Visual_Abstract.png" ]; then
    echo "ERROR: Experiment_Visual_Abstract.png not found in current directory"
    exit 1
fi

echo "✓ All required files found"
echo ""

# Create output directory
mkdir -p "$OUT_DIR"

# ==============================================================================
# CONTAINER SETUP
# ==============================================================================
IMG_PATH="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"

if [[ ! -f "$IMG_PATH" ]]; then
    echo "Pulling Bioconductor container..."
    singularity pull "$IMG_PATH" docker://go2432/bioconductor:latest
fi

echo "✓ Container ready: $IMG_PATH"
echo ""

# ==============================================================================
# RUN PUBLICATION FIGURE GENERATION
# ==============================================================================
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║      GENERATING PUBLICATION FIGURE (Panels A-I)                ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "Configuration:"
echo "  VST File:         $VST_FILE"
echo "  Metadata File:    $METADATA_FILE"
echo "  Results Dir:      $RESULTS_DIR"
echo "  GMT Dir:          $GMT_DIR"
echo "  Combined GMT:     $COMBINED_GMT"
echo "  STRING Dir:       $STRING_DIR"
echo "  Output Dir:       $OUT_DIR"
echo "  Target Contrast:  $TARGET_CONTRAST"
echo ""

singularity exec --bind "$PWD:/data" --pwd /data "$IMG_PATH" \
    Rscript /data/create_publication_figure.R \
    "/data/$VST_FILE" \
    "/data/$RESULTS_DIR" \
    "/data/$GMT_DIR" \
    "/data/$STRING_DIR" \
    "/data/$OUT_DIR" \
    "$TARGET_CONTRAST" \
    "/data/$METADATA_FILE"

exit_code=$?

if [ $exit_code -ne 0 ]; then
    echo ""
    echo "ERROR: Publication figure generation failed (exit code: $exit_code)"
    exit $exit_code
fi

# ==============================================================================
# VERIFY OUTPUTS
# ==============================================================================
echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║      VERIFYING OUTPUTS                                         ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

required_outputs=(
    "$OUT_DIR/Publication_Figure_9Panel.png"
    "$OUT_DIR/Publication_Figure_9Panel.pdf"
    "$OUT_DIR/Figure_Captions.txt"
)

all_present=true
for output in "${required_outputs[@]}"; do
    if [ -f "$output" ]; then
        size=$(du -h "$output" | cut -f1)
        echo "✓ $output ($size)"
    else
        echo "✗ MISSING: $output"
        all_present=false
    fi
done

if [ "$all_present" = false ]; then
    echo ""
    echo "WARNING: Some expected outputs are missing"
    exit 1
fi

# ==============================================================================
# SUMMARY
# ==============================================================================
echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║      PUBLICATION FIGURE GENERATION COMPLETE                    ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "📊 GENERATED OUTPUTS:"
echo ""
echo "Main Figure:"
echo "  • PNG (high-res): $OUT_DIR/Publication_Figure_9Panel.png"
echo "  • PDF (vector):   $OUT_DIR/Publication_Figure_9Panel.pdf"
echo "  • Captions:       $OUT_DIR/Figure_Captions.txt"
echo ""
echo "Panel Layout (3×3 grid):"
echo "  A. Experimental Design (Gemini visual abstract)"
echo "  B. Global Structure (PCA biplot + scree plot)"
echo "  C. Subtype Trajectories (limma pairwise with arrayWeights)"
echo "  D. Semantic Pathway Clustering (tree plot)"
echo "  E. Protein-Protein Interaction Network"
echo "  F. Polypharmacology Network"
echo "  G. Drug-Pathway Gene Overlap Heatmap"
echo "  H. Drug Candidates 2D Plot (BBB × NES)"
echo "  I. Top 15 Drug Candidates Table"
echo ""
echo "📝 PANEL DESCRIPTIONS:"
echo ""
cat "$OUT_DIR/Figure_Captions.txt"
echo ""
echo "════════════════════════════════════════════════════════════════"
echo "Statistical Methods (Panel C):"
echo "  • Limma linear modeling with arrayWeights"
echo "  • All pairwise contrasts tested"
echo "  • Significance: * p<0.05, ** p<0.01, *** p<0.001"
echo "════════════════════════════════════════════════════════════════"
echo "                    ANALYSIS COMPLETE"
echo "════════════════════════════════════════════════════════════════"
