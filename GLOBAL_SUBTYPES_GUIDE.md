# Global Subtypes Analysis v16.1 - User Guide

## 🎯 Purpose
Analyze GBM subtype evolution after Litt therapy by detecting **trajectory patterns** across Culture → Primary → Recurrent progression.

## 🆕 What's New in v16.1

### 1. **Automatic Dual-Weighting Comparison**
- Runs **TWO analyses automatically** (no flags needed):
  - **Weighted**: Uses arrayWeights (standard, downweights outliers)
  - **Unweighted**: All samples equal (unbiased, preserves subtle switches)
- **Why this matters**: arrayWeights can mistake biological switches for "noise"
- **How to interpret**: 
  - Significant in BOTH → Robust finding ✓
  - Significant ONLY in unweighted → arrayWeights masked it (likely real!)

### 2. **Dual Plasticity Testing (ANOVA + Kruskal-Wallis)**
- **Runs BOTH tests automatically** and displays results side-by-side
- **Why both?**
  - ANOVA: More powerful when data is normal
  - Kruskal-Wallis: More robust for small N or non-normal data
- **When they agree**: Highest confidence in results
- **When they disagree**: Check normality test to decide which to trust

### 3. **Trajectory-Aware Statistical Testing**
Classic pairwise tests (ANOVA, limma) ignore your experimental design's **ordered nature**.

v16.1 adds THREE trajectory tests:

#### A. Jonckheere-Terpstra Test
- **What it does**: Tests for monotonic trends (steady increase/decrease)
- **Why use it**: More powerful than ANOVA for ordered groups
- **Example**: If Neftel_AC goes 1.0 → 1.5 → 2.0, JT says "upward trajectory" (not just "different")

#### B. Polynomial Contrasts (Linear)
- **What it does**: Tests for steady progression
- **Interpretation**: 
  - Linear P < 0.05 = Signature steadily increases/decreases
  - Example: Mesenchymal drift across evolution

#### C. Polynomial Contrasts (Quadratic)
- **What it does**: Tests for "spike/dip" patterns
- **Interpretation**:
  - Quadratic P < 0.05 = Signature peaks or valleys in middle stage
  - Example: Therapy shock response (high in Primary, low in Culture & Recurrent)

### 3. **Pattern Classification**
Each signature is automatically classified:
- **Linear (monotonic)**: JT P < 0.05, linear P < 0.05, quadratic ns
- **Quadratic (spike/dip)**: Quadratic P < 0.05
- **Weak trend**: JT P < 0.10
- **No trend**: None of the above

### 4. **LLM-Ready Interpretation Summary**
Creates `_llm_summary.txt` with:
- Comparison of weighted vs unweighted results
- List of significant trajectory patterns
- Interpretation guide for biological context
- **Usage**: Copy file contents into ChatGPT/Claude/Gemini

## 📊 Usage

### Basic Command
```bash
Rscript run_global_subtypes_v16_1.R \
  -i vst_expression.tsv \
  -m metadata.csv \
  -o results/litt_evolution
```

### Required Input Files

**1. Expression Matrix** (`-i`)
- Format: TSV with genes as rows, samples as columns
- Can be Ensembl IDs (auto-converted) or gene symbols
- Example:
```
GeneID    Sample1    Sample2    Sample3
ENSG00001    12.5      13.2      11.8
ENSG00002    8.3       9.1       7.9
```

**2. Metadata** (`-m`)
- Format: CSV with rownames = sample IDs
- **Required column**: `Classification`
- Classification must have ordered levels: `Culture_U2`, `Primary_U2`, `Recurrent_U2`
- Example:
```
Sample,Classification
Sample1,Culture_U2
Sample2,Primary_U2
Sample3,Recurrent_U2
```

### Optional Parameters
- `--n_top_genes 500`: Number of variable genes for PCA (default: 500)
- `--scoring both`: Method for signature scoring (zscore/gsva/both)
- `--signatures custom.gmt`: Use custom signature set instead of built-in

## 📁 Output Files

### MultiQC HTML Report (MAIN OUTPUT)
**`analysis_summary_mqc.html`** ⭐ **VIEW THIS FIRST**
   - Complete analysis results in one place
   - **Embedded LLM summary** (no need for separate file!)
   - Plasticity test comparison table (ANOVA vs Kruskal-Wallis)
   - Trajectory pattern summary
   - All statistical results with explanations
   
### Key CSV Results
1. **`_Trajectory_Tests.csv`**
   - Jonckheere-Terpstra P-values
   - Linear/Quadratic P-values
   - Pattern classification
   
2. **`_Plasticity_Tests_Comparison.csv`**
   - Shapiro-Wilk normality test
   - ANOVA results
   - Kruskal-Wallis results
   - Side-by-side comparison
   
3. **`_llm_summary.txt`** (also embedded in HTML)
   - Standalone LLM interpretation prompt
   - Copy into ChatGPT/Claude/Gemini
   - Compares weighted vs unweighted findings

### Standard Outputs
- `_Global_Structure_mqc.png`: PCA with biplot
- `_Heatmap_mqc.png`: Signature heatmap
- `_Plasticity_mqc.png`: Shannon entropy by group
- `_Scores.csv`: Signature scores per sample
- `_Metadata.csv`: Metadata with plasticity scores

## 🔍 Interpreting Results

### Step 1: Open the MultiQC HTML Report
**File**: `analysis_summary_mqc.html`

This report contains:
- ✅ Global structure (PERMANOVA, PCA)
- ✅ Plasticity test comparison table (ANOVA vs KW)
- ✅ Trajectory pattern summary
- ✅ Top significant signatures
- ✅ **Embedded LLM interpretation summary** (scroll down!)

### Step 2: Check Plasticity Tests
Look at the "Plasticity Analysis (Dual Testing)" table:
- **Shapiro-Wilk**: Is data normal? (P < 0.05 = non-normal)
- **ANOVA**: Parametric test result
- **Kruskal-Wallis**: Non-parametric test result

**Interpretation**:
- Both significant (P < 0.05) → High confidence ✓
- Both non-significant → No plasticity differences
- Disagree → Use test recommended at bottom of table

### Step 3: Identify True Trajectory Signatures
Look at the "Top Trajectory Signatures" table in HTML report, OR:
```r
library(dplyr)
traj <- read.csv("results/litt_evolution_Trajectory_Tests.csv")

# Find strongest trajectory signals
traj %>%
  filter(JT_P < 0.05) %>%
  arrange(JT_P) %>%
  select(Signature, JT_Direction, JT_P, Pattern)
```

### Step 4: Compare Weighting Methods
Scroll to the "LLM Interpretation Summary" section in the HTML report (or open `_llm_summary.txt`):
- Signatures significant in BOTH → High confidence
- Signatures ONLY in unweighted → arrayWeights masked (still likely real)

### Step 5: Get AI Interpretation
**Option A (Easiest)**: Copy the "LLM Interpretation Summary" section from the HTML report

**Option B**: Copy the entire `_llm_summary.txt` file

Paste into ChatGPT/Claude/Gemini with:
```
Interpret these GBM subtype evolution results in the context of Litt therapy.
Focus on:
1. Which subtypes show clear evolutionary trajectories?
2. What do linear vs quadratic patterns suggest biologically?
3. Clinical implications for recurrence mechanisms?
4. Suggested follow-up experiments?
```

## 🎓 Understanding the Statistics

### Why Jonckheere-Terpstra?
Standard tests (ANOVA, Kruskal-Wallis) ask: "Are any groups different?"
JT asks: "Is there an ordered trend?"

Example:
- Values: Culture=1.0, Primary=1.2, Recurrent=2.0
- ANOVA: P=0.03 (groups differ)
- JT: P=0.008 (stronger, detects the ORDER matters)

### Why Polynomial Contrasts?
They test **HOW** a signature changes:
- **Linear**: Steady progression (e.g., EMT: 1.0 → 1.5 → 2.0)
- **Quadratic**: Transient response (e.g., stress: 1.0 → 3.0 → 1.2)

### arrayWeights Issue
In small datasets (N < 10/group), biological switches can look like outliers:
- Culture: 1.0, 1.1, 1.0
- Primary: 1.2, 1.3, 1.1
- Recurrent: 2.5, 2.6, 2.4  ← arrayWeights thinks these are "bad samples"

Solution: v16.1 runs BOTH analyses so you can compare

## 📚 Built-in Signatures
- **Verhaak (4)**: Classical, Mesenchymal, Proneural, Neural
- **Neftel (4)**: AC, OPC, NPC, MES
- **Garofano (3)**: MTC, GPM, NEU

Total: 11 GBM subtype signatures

## 🐛 Troubleshooting

### Error: "No genes found for signature X"
- Check gene IDs match (Ensembl vs symbols)
- Verify expression matrix has gene names as rownames

### Warning: "arrayWeights convergence failed"
- Normal for very small datasets
- Results still valid (uses equal weights as fallback)

### All JT P-values > 0.05
- Stages may not have clear trajectory
- Check PCA plot: do stages separate?
- Try `--n_top_genes 1000` for more features

## 📧 Support
See full code comments in `run_global_subtypes_v16_1.R` for implementation details.

---
**Version**: 16.1 (Dual-Weighting + Trajectory Edition)  
**Date**: January 2026  
**Optimized for**: Litt therapy evolution studies
