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

# ==============================================================================
# SETUP
# ==============================================================================
export CONDA_PREFIX="${HOME}/mambaforge/envs/nextflow"
export PATH="${CONDA_PREFIX}/bin:$PATH"
unset JAVA_HOME
nextflow -v

# Cache directories
export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
mkdir -p $XDG_RUNTIME_DIR $NXF_SINGULARITY_CACHEDIR

# Safety - clean environment
export NXF_SINGULARITY_HOME_MOUNT=true
unset LD_LIBRARY_PATH
unset PYTHONPATH
unset R_LIBS
unset R_LIBS_USER
unset R_LIBS_SITE

# ==============================================================================
# PATH CONFIGURATION
# ==============================================================================
GMT_DIR="ANALYSIS/refs/pathways/human"
STRING_DIR="ANALYSIS/refs/pathways/human/human_string"
RESULTS_DIR="ANALYSIS/results_evolution"
VST_FILE="${RESULTS_DIR}/tables/processed_abundance/all.vst.tsv"

# Output directory for MultiQC-ready files
OUT_DIR="ANALYSIS/results_visualization/Ultimate_Report"
mkdir -p "$OUT_DIR"

# Singularity containers
IMG_PATH="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"
MULTIQC_CONTAINER="docker://multiqc/multiqc:v1.33"

if [[ ! -f "$IMG_PATH" ]]; then
    echo "Pulling Bioconductor container..."
    singularity pull "$IMG_PATH" docker://go2432/bioconductor:latest
fi

# ==============================================================================
# STEP 1: ULTIMATE PATHWAY/DRUG/PPI ANALYSIS
# ==============================================================================
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║  STEP 1: Ultimate Pathway Enrichment & Drug Discovery         ║"
echo "╚════════════════════════════════════════════════════════════════╝"

singularity exec --bind $PWD:/data --pwd /data "$IMG_PATH" \
    Rscript run_pathways_drugs_v6_SUPREME.R \
    "$VST_FILE" \
    "$RESULTS_DIR" \
    "$GMT_DIR" \
    "$STRING_DIR" \
    "$OUT_DIR/Analysis"

if [ $? -ne 0 ]; then
    echo "ERROR: Pathway analysis failed"
    exit 1
fi

echo "✓ Pathway analysis complete"
echo ""

# ==============================================================================
# STEP 2: GENERATE MULTIQC REPORT
# ==============================================================================
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║  STEP 2: MultiQC Aggregation                                  ║"
echo "╚════════════════════════════════════════════════════════════════╝"

# Create MultiQC configuration
cat << 'CONFIG' > mqc_config.yaml
# MultiQC Configuration for Ultimate Analysis Report
title: "U251 Transcriptomic Evolution - Ultimate Analysis"
subtitle: "Pathway Enrichment, Drug Discovery, and PPI Networks"
intro_text: "Comprehensive analysis of GBM transcriptional evolution with therapeutic target identification"

# Disable version detection to reduce clutter
disable_version_detection: true

# Hide sections we don't need
remove_sections:
  - software_versions

# Custom content ordering
custom_content:
  order:
    - analysis_narrative
    - pathway_enrichment
    - drug_discovery
    - ppi_networks

# Search patterns for custom content
sp:
  analysis_narrative:
    fn: "*Narrative_mqc.html"
  pathway_enrichment:
    fn: "*GSEA*_mqc.png"
  drug_discovery:
    fn: "*Drug*_mqc.png"
  ppi_networks:
    fn: "*PPI*_mqc.png"

# Report appearance
report_header_info:
  - Contact E-mail: 'analysis@institution.edu'
  - Application Type: 'RNA-seq pathway analysis'
  - Project Type: 'GBM therapeutic target discovery'

# Table configuration
table_columns_visible:
  general_stats: True
CONFIG

# Run MultiQC
singularity exec --bind $PWD:/data --pwd /data $MULTIQC_CONTAINER multiqc \
    /data/$OUT_DIR \
    --force \
    --config /data/mqc_config.yaml \
    --title "U251 LITT Therapy Evolution - Ultimate Analysis" \
    --filename "U251_Ultimate_Report.html" \
    --outdir "/data/$OUT_DIR" \
    --comment "Complete pathway enrichment analysis with drug discovery and PPI network analysis"

if [ $? -ne 0 ]; then
    echo "WARNING: MultiQC failed, but analysis results are still available"
fi

echo "✓ MultiQC report complete"
echo ""

# ==============================================================================
# STEP 3: GENERATE SUMMARY
# ==============================================================================
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║  ANALYSIS COMPLETE - File Summary                             ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "📊 PRIMARY OUTPUTS:"
echo "   • MultiQC Report:        $OUT_DIR/U251_Ultimate_Report.html"
echo "   • Narrative (HTML):      $OUT_DIR/Analysis_Narrative_mqc.html"
echo "   • LLM Prompt (TXT):      $OUT_DIR/LLM_Analysis_Prompt.txt"
echo ""
echo "🔬 ANALYSIS FEATURES:"
echo "   ✓ Complete GSEA pathway enrichment"
echo "   ✓ Drug discovery predictions (DSigDB)"
echo "   ✓ PPI network analysis with hub genes"
echo "   ✓ Comprehensive plot interpretation guides"
echo "   ✓ AI-ready structured summaries"
echo ""
echo "📈 PLOT TYPES GENERATED:"
echo "   • Dotplots (pathway ranking)"
echo "   • Enrichment maps (pathway networks)"
echo "   • Running scores (GSEA validation)"
echo "   • Gene-pathway networks (driver identification)"
echo "   • Ridgeplots (expression distributions)"
echo "   • PPI networks (protein interactions)"
echo ""
echo "🤖 AI/LLM INTEGRATION:"
echo "   • Self-contained HTML with embedded interpretation"
echo "   • Structured TXT prompt for ChatGPT/Claude/Gemini"
echo "   • Complete statistical documentation"
echo "   • Biological context for all analyses"
echo ""

# Count generated files
n_plots=$(find "$OUT_DIR" -name "*_mqc.png" | wc -l)
n_csv=$(find "$OUT_DIR" -name "*.csv" | wc -l)

echo "📁 FILE COUNTS:"
echo "   • Plots generated: $n_plots"
echo "   • Data tables: $n_csv"
echo ""
echo "════════════════════════════════════════════════════════════════"
echo "                    PIPELINE COMPLETE"
echo "════════════════════════════════════════════════════════════════"
