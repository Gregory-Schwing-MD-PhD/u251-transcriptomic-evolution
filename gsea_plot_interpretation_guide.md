# 🔬 GSEA Plot Interpretation Guide
## Complete Guide to Understanding Your Pathway Enrichment Results

---

## Table of Contents
1. [Overview](#overview)
2. [Dotplot](#dotplot)
3. [Enrichment Map (Emap)](#enrichment-map)
4. [Running Enrichment Score](#running-enrichment-score)
5. [Gene-Pathway Network (Cnetplot)](#gene-pathway-network)
6. [Ridgeplot](#ridgeplot)
7. [Drug Discovery Plots](#drug-discovery)
8. [PPI Network](#ppi-network)
9. [Statistical Concepts](#statistical-concepts)
10. [Practical Workflow](#practical-workflow)

---

## Overview

### What is GSEA?
**Gene Set Enrichment Analysis (GSEA)** determines whether a predefined set of genes (pathway, biological process, drug signature) shows statistically significant, coordinated differences between two biological states.

### Key Metrics

#### NES (Normalized Enrichment Score)
- **What it is:** Measure of pathway enrichment accounting for gene set size
- **Range:** Typically -3 to +3 (can be higher)
- **Interpretation:**
  - **Positive NES (Red):** Pathway is ACTIVATED/UPREGULATED in your test condition
  - **Negative NES (Blue):** Pathway is SUPPRESSED/DOWNREGULATED in your test condition
  - **|NES| > 1.5:** Moderate enrichment
  - **|NES| > 2.0:** Strong enrichment
  - **|NES| > 2.5:** Very strong enrichment

#### FDR (False Discovery Rate)
- **What it is:** Probability that this result is a false positive (multiple testing corrected)
- **Interpretation:**
  - **FDR < 0.05:** Statistically significant ✓
  - **FDR 0.05-0.25:** Suggestive/exploratory (use with caution)
  - **FDR > 0.25:** Not significant

#### p.adjust (Adjusted P-value)
- Same as FDR - these terms are used interchangeably
- Benjamini-Hochberg correction for multiple testing

---

## Dotplot

### What It Shows
**Top enriched pathways ranked by statistical significance**, separated into activated and suppressed panels.

### Visual Elements

#### X-axis: GeneRatio
- **Definition:** Proportion of genes in the pathway that are in your significant DE gene list
- **Example:** GeneRatio = 0.25 means 25% of pathway genes are differentially expressed
- **Interpretation:**
  - **Higher GeneRatio (right):** More pathway coverage, stronger biological signal
  - **Lower GeneRatio (left):** Fewer genes, may be driven by specific components

#### Dot Size: Count
- **Definition:** Absolute number of your DE genes in this pathway
- **Interpretation:**
  - **Large dots:** Many genes (more robust finding)
  - **Small dots:** Few genes (may be specific subcomponent)

#### Dot Color: p.adjust (FDR)
- **Scale:** Red (significant) → Blue/Purple (less significant)
- **Interpretation:**
  - **Deep red (p.adjust < 0.01):** High confidence
  - **Yellow/Orange (p.adjust 0.01-0.05):** Significant
  - **Blue/Purple (p.adjust > 0.05):** Not significant

#### Panels: Activated vs Suppressed
- **Left panel:** Pathways with positive NES (upregulated)
- **Right panel:** Pathways with negative NES (downregulated)

### How to Read

**Example:**
```
Pathway: "Oxidative Phosphorylation"
GeneRatio: 0.35
Count: 28 genes
p.adjust: 0.001
Panel: Activated (left)
```

**Interpretation:**
> "The Oxidative Phosphorylation pathway is significantly ACTIVATED. 35% of genes in this pathway (28 total genes) are upregulated in my test condition. This is a high-confidence finding (FDR < 0.001)."

### What to Look For
1. **Top pathways:** Those at the top are most significant
2. **Large red dots:** High confidence findings with many genes
3. **Biological coherence:** Do related pathways cluster together?
4. **Balance:** Are more pathways activated or suppressed?

---

## Enrichment Map

### What It Shows
**Network visualization of pathway relationships based on gene overlap**. Pathways that share many genes are connected.

### Visual Elements

#### Nodes (Circles)
- **Each node = one enriched pathway**
- **Node size:** Number of genes in pathway (larger = more genes)
- **Node color:**
  - **Red/Orange:** Positive NES (activated)
  - **Blue/Purple:** Negative NES (suppressed)
  - **Color intensity:** Magnitude of NES (darker = stronger)

#### Edges (Connections)
- **Definition:** Lines connecting pathways that share genes
- **Edge thickness:** Degree of gene overlap
  - **Thick lines:** Many shared genes
  - **Thin lines:** Few shared genes
- **No connection:** Pathways share few/no genes

#### Network Clusters
- **Definition:** Groups of tightly connected nodes
- **Interpretation:** Functionally related processes

### How to Read

**Example:**
```
Cluster of 5 connected red nodes:
- "Cell Cycle"
- "DNA Replication"
- "Mitotic Spindle"
- "E2F Targets"
- "G2M Checkpoint"
All connected by thick edges
```

**Interpretation:**
> "There is a coordinated ACTIVATION of the proliferation program. Multiple cell cycle pathways are upregulated and share many genes, suggesting a coherent biological response rather than isolated changes."

### What to Look For

1. **Hub pathways:** Nodes with many connections = central processes
2. **Functional modules:** Clusters represent biological themes
3. **Isolated pathways:** May represent unique processes
4. **Color patterns:** Are clusters uniform in color (coordinated) or mixed?
5. **Network density:** Highly connected = many overlapping pathways

### Biological Interpretation Guide

**Dense red cluster:**
> "Coordinated activation of related processes - strong biological program"

**Isolated blue node:**
> "Specific suppression of unique pathway - may be regulatory mechanism"

**Mixed color cluster:**
> "Complex regulation - some components activated, others suppressed within same functional group"

---

## Running Enrichment Score

### What It Shows
**Visual proof of pathway enrichment** showing how genes in a pathway are distributed across your ranked gene list.

### Three Panel Structure

#### Top Panel: Enrichment Score (ES) Curve
- **Y-axis:** Running enrichment score
- **X-axis:** Gene rank (left = most upregulated, right = most downregulated)
- **The curve:**
  - Increases when pathway genes are encountered
  - Decreases when non-pathway genes are encountered
  
**Peak interpretation:**
- **Peak on LEFT:** Pathway genes are enriched among upregulated genes → Pathway ACTIVATED
- **Peak on RIGHT:** Pathway genes are enriched among downregulated genes → Pathway SUPPRESSED
- **Flat curve:** No enrichment - genes evenly distributed

**Peak magnitude (height):**
- **High peak (|ES| > 0.5):** Strong enrichment
- **Low peak (|ES| < 0.3):** Weak enrichment

#### Middle Panel: Gene Hits (Black Bars)
- **Each black line:** Position where a pathway gene appears in ranked list
- **Interpretation:**
  - **Clustered on LEFT:** Pathway upregulated
  - **Clustered on RIGHT:** Pathway downregulated
  - **Evenly spread:** No enrichment (false positive)

#### Bottom Panel: Ranking Metric
- **Usually log2 Fold Change or Wald statistic**
- **Shows actual gene-level statistics**
- **Context for understanding enrichment**

### How to Read

**Example: "Good" Enrichment**
```
Top panel: Peak at position 2000 (left side), ES = 0.65
Middle panel: Most black bars clustered between rank 0-5000
Bottom panel: High positive fold changes on left
```

**Interpretation:**
> "STRONG ACTIVATION. Pathway genes are concentrated among the most upregulated genes, not scattered throughout. This is a robust finding."

**Example: "Bad" Enrichment (False Positive)**
```
Top panel: Peak at position 8000, ES = 0.42
Middle panel: Black bars evenly spread across all ranks
Bottom panel: Mixed fold changes
```

**Interpretation:**
> "WEAK/QUESTIONABLE ENRICHMENT. Even though GSEA reports significance, genes are not well-clustered. This may be a false positive driven by a few outlier genes."

### What to Look For

**✅ Signs of TRUE enrichment:**
- Clear peak in ES curve
- Gene hits clustered at one end
- Peak aligns with gene hit concentration
- Many genes contributing (dense black bars)

**❌ Signs of FALSE POSITIVE:**
- Weak/broad peak
- Gene hits evenly distributed
- Few gene hits
- Peak doesn't align with hit concentration

---

## Gene-Pathway Network

### What It Shows
**Which specific genes are driving pathway enrichment** and how genes connect multiple pathways.

### Visual Elements

#### Large Nodes: Pathways
- Outer circle nodes represent enriched pathways
- Size indicates number of genes

#### Small Nodes: Genes
- Inner nodes represent individual genes
- **Color scale:**
  - **Deep red:** Highly upregulated (large positive log2FC)
  - **Light red:** Moderately upregulated
  - **Light blue:** Moderately downregulated
  - **Deep blue:** Highly downregulated (large negative log2FC)

#### Edges: Gene-Pathway Membership
- Lines connect genes to pathways they belong to
- A gene connected to multiple pathways is a "hub gene"

### How to Read

**Example:**
```
Gene: "TP53"
Color: Deep red
Connected to: "Apoptosis", "DNA Damage Response", "Cell Cycle"
```

**Interpretation:**
> "TP53 is highly upregulated (deep red) and participates in multiple cancer-related pathways. It's a key regulatory node driving enrichment across several processes."

### What to Look For

1. **Hub genes:** Connected to many pathways
   - Often transcription factors, signaling proteins
   - Potential therapeutic targets
   - Master regulators

2. **Extreme colors:** Genes with largest fold changes
   - Deep red/blue = strongest expression changes
   - These are your top hits

3. **Pathway-specific genes:** Connected to only one pathway
   - Specialized function
   - May represent specific mechanisms

4. **Gene-pathway coherence:**
   - Do genes in same pathway have similar colors?
   - Uniform color = coordinated regulation
   - Mixed colors = complex regulation

### Biological Use Cases

**Finding drug targets:**
> Look for hub genes with extreme colors (deep red/blue) that connect many cancer-related pathways

**Understanding mechanism:**
> Trace from pathway → genes to identify which specific genes are driving enrichment

**Validation priorities:**
> Hub genes with extreme fold changes are top candidates for qPCR validation

---

## Ridgeplot

### What It Shows
**Distribution of gene expression changes within each enriched pathway**, revealing whether regulation is uniform or mixed.

### Visual Elements

#### Each Ridge (Row)
- One enriched pathway
- Shows distribution of log2FC values for genes in that pathway

#### X-axis: log2 Fold Change
- Negative values (left): Downregulation
- Positive values (right): Upregulation
- Zero: No change

#### Ridge Shape & Position

**Shifted RIGHT (positive):**
- Peak is in positive territory
- Most pathway genes are upregulated
- Confirms pathway activation

**Shifted LEFT (negative):**
- Peak is in negative territory
- Most pathway genes are downregulated
- Confirms pathway suppression

**Narrow, tall peak:**
- Genes have similar fold changes
- Coordinated, uniform regulation
- High confidence

**Wide, flat peak:**
- Genes have variable fold changes
- Mixed regulation
- May indicate pathway complexity

**Bimodal (two peaks):**
- Some genes up, some down within same pathway
- Complex regulation
- Pathway has opposing components

### How to Read

**Example 1: "Clean" Activation**
```
Pathway: "Interferon Response"
Ridge: Narrow peak at +2.5 log2FC
Width: Very tight distribution
```

**Interpretation:**
> "STRONG UNIFORM ACTIVATION. Nearly all interferon response genes are upregulated by ~2.5-fold. This is a coordinated biological response."

**Example 2: "Complex" Regulation**
```
Pathway: "TGF-beta Signaling"
Ridge: Bimodal - one peak at -1.5, another at +2.0
Width: Very wide
```

**Interpretation:**
> "COMPLEX REGULATION. Some TGF-beta pathway components are suppressed while others are activated. This suggests selective pathway modulation rather than complete activation/suppression."

### What to Look For

**✅ High confidence findings:**
- Narrow peaks (consistent regulation)
- Peak far from zero (strong effect)
- Unimodal distribution (one clear peak)

**⚠️ Complex regulation:**
- Wide distributions (variable effects)
- Bimodal peaks (opposing regulation)
- Peak near zero (weak effects)

---

## Drug Discovery

### Special Interpretation for Drug Signatures

**KEY CONCEPT:** For drug discovery, we want drugs whose signatures are **OPPOSITE** to disease.

#### NES Interpretation for Drugs

**Negative NES (Blue) = THERAPEUTIC POTENTIAL ✓**
- Drug signature opposes disease signature
- Drug downregulates what disease upregulates
- Drug upregulates what disease downregulates
- **More negative NES = stronger therapeutic potential**

**Positive NES (Red) = AVOID ⚠️**
- Drug signature mimics disease signature
- Drug may worsen disease
- Not a therapeutic candidate

### How to Prioritize Drug Candidates

**Tier 1 (Highest Priority):**
```
NES < -2.0
FDR < 0.01
FDA-approved status
```

**Tier 2 (Moderate Priority):**
```
NES < -1.5
FDR < 0.05
Clinical trial status
```

**Tier 3 (Exploratory):**
```
NES < -1.0
FDR < 0.25
Preclinical compounds
```

### Drug Plot Interpretation

**Example:**
```
Drug: "DOXORUBICIN"
NES: -2.3
FDR: 0.003
Core enrichment: 45 genes
```

**Interpretation:**
> "Doxorubicin shows STRONG therapeutic potential (NES = -2.3, p < 0.003). Its gene expression signature is strongly opposite to the disease signature, suggesting it would reverse pathological gene expression. This is a high-priority candidate for validation."

### Important Caveats

⚠️ **This is computational prediction only**
- Requires experimental validation
- Does NOT account for:
  - Drug toxicity
  - Drug delivery/bioavailability
  - Pharmacokinetics
  - Off-target effects
  - Clinical feasibility

---

## PPI Network

### What It Shows
**Protein-protein interaction network** of your significant genes based on STRING database.

### Visual Elements

#### Hub Genes (Red, Large)
- Proteins with many interactions (high degree)
- Often transcription factors, kinases, signaling proteins
- **Interpretation:**
  - Master regulators
  - Central to cellular processes
  - Disrupting these affects many downstream pathways
  - Prime therapeutic targets

#### Regular Nodes (Blue, Small)
- Less connected proteins
- May be specialized effectors

#### Network Clusters
- Groups of tightly connected proteins
- Often represent:
  - Protein complexes (e.g., ribosome, proteasome)
  - Functional modules (e.g., DNA repair machinery)
  - Signaling cascades

### How to Read

**Example:**
```
Hub Gene: "MYC"
Connections: 25 proteins
Cluster: Central node connecting cell cycle and metabolic clusters
```

**Interpretation:**
> "MYC is a master regulator (hub) connecting proliferation and metabolism. It's differentially expressed and central to the network, making it a high-priority target for validation and potential therapeutic intervention."

### What to Look For

1. **Hub genes:**
   - Top candidates for validation
   - Potential biomarkers
   - Drug target opportunities

2. **Functional modules:**
   - Clusters of connected genes
   - Represent coordinated processes
   - May indicate:
     - Protein complexes to target
     - Pathway dependencies
     - Synthetic lethal opportunities

3. **Network topology:**
   - **Dense, interconnected:** Many regulatory relationships
   - **Sparse, modular:** Independent processes
   - **Star pattern:** One master regulator with many targets

### Therapeutic Implications

**Hub genes = therapeutic targets:**
```
Rationale:
1. Affect many pathways (broad impact)
2. Often druggable (kinases, TFs with small molecule inhibitors)
3. May have synthetic lethal partners
4. Potential biomarkers (centrality suggests importance)
```

**Follow-up experiments:**
1. **Validate expression:** qPCR, Western blot
2. **Functional testing:** siRNA knockdown, CRISPR knockout
3. **Drug screening:** Existing inhibitors in DrugBank
4. **Mechanistic studies:** ChIP-seq (if TF), phosphoproteomics (if kinase)

---

## Statistical Concepts

### Multiple Testing Correction

**Why it's needed:**
- Testing thousands of pathways simultaneously
- Some will appear significant by chance (false positives)
- Need to control false discovery rate

**Benjamini-Hochberg (FDR) correction:**
- Standard approach in GSEA
- Controls expected proportion of false discoveries
- More powerful than Bonferroni (less conservative)

### Permutation Testing

**How GSEA computes p-values:**
1. Shuffle gene labels randomly (1000 times)
2. Recompute enrichment score each time
3. Compare real ES to null distribution
4. p-value = fraction of permutations with ES ≥ observed

**Advantages:**
- Non-parametric (no distribution assumptions)
- Accounts for gene-gene correlation
- Robust to outliers

### Gene Set Size

**Why it matters:**
- Very small sets (<15 genes): Unstable statistics
- Very large sets (>500 genes): Loss of specificity
- Optimal: 15-500 genes

**Default thresholds in this pipeline:**
```R
minGSSize = 15
maxGSSize = 500
```

---

## Practical Workflow

### Step 1: Start with Dotplot
**Goal:** Get overview of top pathways

**Questions:**
- Which pathways are most significant?
- Are more pathways activated or suppressed?
- Do I see biological coherence (related pathways together)?

**Action:**
- Identify top 5-10 pathways in each direction
- Screenshot for presentation

---

### Step 2: Validate with Running Score
**Goal:** Confirm enrichment is real, not statistical artifact

**Questions:**
- Are gene hits well-clustered?
- Is the peak clear and strong?
- Do many genes contribute?

**Action:**
- Check running score for top 3-5 pathways
- Discard any with poor clustering
- Keep only those with clear peaks

---

### Step 3: Explore Relationships with Emap
**Goal:** Understand pathway crosstalk

**Questions:**
- Do my top pathways form clusters?
- Are there hub pathways connecting many others?
- Does the network tell a biological story?

**Action:**
- Identify functional modules (clusters)
- Note hub pathways
- Draw conceptual model of pathway relationships

---

### Step 4: Identify Driver Genes with Network
**Goal:** Find specific genes causing enrichment

**Questions:**
- Which genes are hubs (connect many pathways)?
- Which genes have extreme fold changes?
- Are there druggable targets?

**Action:**
- List top 10 hub genes
- Note genes with |log2FC| > 2
- Cross-reference with DrugBank
- Prioritize for validation

---

### Step 5: Assess Regulation with Ridgeplot
**Goal:** Understand consistency of regulation

**Questions:**
- Is regulation uniform or mixed?
- Are there bimodal patterns suggesting complexity?
- Do distributions match my expectations?

**Action:**
- Note pathways with narrow peaks (uniform)
- Flag bimodal patterns for further investigation
- Use for mechanistic interpretation

---

### Step 6: Drug Discovery (if applicable)
**Goal:** Identify therapeutic candidates

**Questions:**
- Which drugs oppose disease signature?
- Are they FDA-approved or in trials?
- Do they target hub genes or key pathways?

**Action:**
- List drugs with NES < -1.5, FDR < 0.05
- Check FDA status
- Prioritize based on mechanism
- Plan validation experiments

---

### Step 7: PPI Network Analysis
**Goal:** Identify master regulators and targets

**Questions:**
- Which hub genes are central to network?
- Do hubs participate in enriched pathways?
- Are there druggable hubs?

**Action:**
- List top 10 hubs
- Check pathway membership
- Literature review for each
- Plan validation (qPCR, Western, functional assays)

---

## Common Pitfalls & How to Avoid

### Pitfall 1: Over-interpreting FDR 0.05-0.25
**Problem:** These are exploratory, not definitive
**Solution:** Require additional evidence (biological plausibility, validation)

### Pitfall 2: Ignoring Running Score
**Problem:** Trusting GSEA p-value without visual validation
**Solution:** Always check running score for gene clustering

### Pitfall 3: Cherry-picking Pathways
**Problem:** Focusing only on expected results
**Solution:** Report top pathways regardless of expectation

### Pitfall 4: Not Validating Hub Genes
**Problem:** Assuming central genes are important without testing
**Solution:** Experimental validation (qPCR, functional assays)

### Pitfall 5: Over-reliance on Drug Predictions
**Problem:** Computational predictions ≠ clinical efficacy
**Solution:** Treat as hypotheses requiring experimental validation

---

## Example Complete Interpretation

### Scenario: GBM Recurrence vs Primary

#### Dotplot Findings:
```
Top Activated:
1. Epithelial-Mesenchymal Transition (NES=2.8, FDR=0.001)
2. Hypoxia (NES=2.3, FDR=0.003)
3. Angiogenesis (NES=2.1, FDR=0.008)

Top Suppressed:
1. Oxidative Phosphorylation (NES=-2.5, FDR=0.001)
2. DNA Repair (NES=-1.9, FDR=0.015)
```

#### Running Score Validation:
- ✅ EMT: Excellent clustering, peak at rank 500
- ✅ Oxidative Phosphorylation: Clear suppression, peak at rank 15000

#### Emap Insights:
- Dense cluster: EMT + Hypoxia + Angiogenesis all connected
- Isolated: DNA Repair (unique suppression)

#### Network Analysis:
```
Top Hubs:
- TWIST1 (25 connections, log2FC=3.2)
- HIF1A (22 connections, log2FC=2.8)
- VIM (20 connections, log2FC=3.5)
```

#### Drug Candidates:
```
1. Bevacizumab (NES=-2.1, FDR=0.007) - anti-angiogenic
2. Metformin (NES=-1.8, FDR=0.012) - metabolic
```

### Integrated Biological Narrative:

> "GBM recurrence shows coordinated activation of the invasive mesenchymal program (EMT, NES=2.8), driven by hub transcription factors TWIST1 and HIF1A responding to hypoxic stress. This is accompanied by suppression of oxidative metabolism, consistent with Warburg effect and aggressive phenotype. The tight clustering of EMT/hypoxia/angiogenesis pathways suggests a coherent adaptive response to therapy.
>
> Therapeutic strategy: Target the hypoxia-EMT axis with bevacizumab (anti-VEGF) combined with metabolic reprogramming via metformin. Hub genes TWIST1 and HIF1A are priority validation targets for mechanistic studies and potential small molecule inhibitor screening."

---

## Quick Reference Card

### Dotplot
**When:** First overview
**Shows:** Top pathways by significance
**Look for:** Large red dots at top

### Emap
**When:** Understanding relationships
**Shows:** Pathway networks
**Look for:** Functional clusters, hub pathways

### Running Score
**When:** Validation
**Shows:** Proof of enrichment
**Look for:** Clear peak, gene clustering

### Network
**When:** Finding drivers
**Shows:** Gene-pathway connections
**Look for:** Hub genes, extreme colors

### Ridgeplot
**When:** Assessing regulation
**Shows:** FC distributions
**Look for:** Narrow peaks, unimodal

### Drug Plots
**When:** Therapeutic discovery
**Shows:** Drug-disease opposition
**Look for:** Negative NES, low FDR

### PPI Network
**When:** Target identification
**Shows:** Protein interactions
**Look for:** Hub genes, clusters

---

## Additional Resources

### Recommended Reading
1. Subramanian et al. (2005) PNAS - Original GSEA paper
2. Korotkevich et al. (2021) - fgsea implementation
3. Yu et al. (2012) - clusterProfiler documentation

### Online Tools
- GSEA website: https://www.gsea-msigdb.org
- clusterProfiler book: https://yulab-smu.top/biomedical-knowledge-mining-book/
- STRING database: https://string-db.org

### Databases Used
- **MSigDB:** Molecular Signatures Database (pathways)
- **DSigDB:** Drug Signatures Database (therapeutic)
- **STRING:** Search Tool for Retrieval of Interacting Genes/Proteins

---

## Contact & Citations

**Pipeline Author:** [Your contact information]

**Key Citations:**
```
GSEA: Subramanian et al. (2005) PNAS 102(43):15545-50
clusterProfiler: Yu et al. (2012) OMICS 16(5):284-287
STRING: Szklarczyk et al. (2021) Nucleic Acids Res 49(D1):D605-612
```

---

**Document Version:** 1.0 (Ultimate Edition)  
**Last Updated:** 2025-01-31  
**Pipeline Version:** run_pathways_drugs_v5_ULTIMATE.R
