# Transcriptomic Evolution of U251 Glioblastoma Cells: From Culture to Orthotopic Growth and Post-therapy Recurrence

![Visual Abstract](./ASSETS/Visual_abstract.png)

## Project Overview

This project addresses a critical gap in glioblastoma (GBM) modeling: the lack of longitudinal RNA-sequencing (RNAseq) datasets tracking tumor evolution through clinically relevant stages. Specifically, we investigate the transcriptomic shifts of U251 cells from **standard in vitro culture**, to **orthotopic growth** in the murine brain (Primary), and finally through **recurrence following Laser Interstitial Thermal Therapy (LITT)**.

While U251 is a widely used model, classic profiling shows that gene expression shifts substantially when implanted orthotopically, underscoring the dominance of the brain microenvironment. This project aims to quantify these adaptations and define the transcriptional consequences of focal thermal ablation.

---

## Experimental Design & Sample Cohorts

This dataset comprises four distinct biological groups representing the trajectory of tumor evolution and relevant controls.

### 1. In Vitro Culture (Baseline)
* **Sample IDs:** `C2B`
* **Description:** U251 human glioblastoma cells maintained in standard 2D adherent culture (log phase growth) prior to implantation.
* **Purpose:** Serves as the transcriptional baseline to identify genes differentially regulated solely by the transition to the brain microenvironment.

### 2. Primary Orthotopic Xenograft (Pre-LITT)
* **Sample IDs:** `IL64B`, `IL67B`, `IL68B`, `IL69B`
* **Model:** Adult female immunodeficient **RNU/RNU rats**.
* **Implantation:** $5 \times 10^5$ U251 cells stereotactically injected into the **striatum** (Coordinates: 3.5 mm right of bregma, depth 3.0 mm).
* **Tumor Status:** Tumors were allowed to establish and grow for approximately 2 weeks (reaching ~4 mm diameter) as confirmed by MRI and dynamic contrast-enhanced (DCE) imaging prior to intervention.
* **Purpose:** Captures the "brain-adapted" signature, highlighting upregulation of ECM, invasion, and vascular programs absent in plastic culture.

### 3. Recurrent Tumor (Post-LITT)
* **Sample IDs:** `IL66B`, `NL70B`, `NL71B`
* **Ablation Method:** Orthotopic tumors were treated with the **Visualase® clinical LITT system** (Medtronic) using a 980 nm diode laser fiber inserted stereotactically into the tumor center.
* **Parameters:** Ablation was performed under real-time MRI guidance (DWI) using rat-adapted settings (1 Volt, 30–40s duration) to achieve coagulative necrosis of the tumor bulk while creating a sublethal thermal penumbra in the peritumoral tissue.
* **Recurrence:** These samples represent the tumor regrowth harvested longitudinally after the thermal ablation procedure, representing the therapy-resistant subclone or stress-adapted state.

### 4. Control Samples
* **Sample IDs:** `N168B`, `N269B`
* **Description:** Non-tumor brain tissue or procedural controls.
* **Purpose:** Provides a negative control background for normalizing tumor-specific expression and identifying non-specific sequencing artifacts.

### Key Reference for Methodology
The orthotopic tumor model and the specific adaptation of the Visualase LITT system for this dataset are detailed in:

> **Adaptation of laser interstitial thermal therapy for tumor ablation under MRI monitoring in a rat orthotopic model of glioblastoma**
> *Nagaraja TN, Bartlett S, Farmer KG, et al.* (2021)
> [**Read Full Text (PMC)**](https://pmc.ncbi.nlm.nih.gov/articles/PMC8893160/)

---

## Background & Scientific Context

### 1. In Vitro vs. Orthotopic Xenografts
Evidence suggests that the brain microenvironment exerts a dominant influence on transcriptional states.
* **Convergent Profiles:** Gene expression profiles of GBM lines (U251, U87) become more similar to each other—and to patient GBM—in orthotopic settings than in vitro.
* **Key Pathways:** Orthotopic tumors upregulate extracellular matrix (ECM) remodeling, cell adhesion, and angiogenesis, whereas in vitro cultures favor classical proliferation genes.
* **Clinical Relevance:** Orthotopic xenografts cluster closely with patient GBM samples, recapitulating hypoxia and invasion signatures that subcutaneous models fail to capture.

**Key Reference:**
* *Influence of in vivo growth on human glioma cell line gene expression* (De Witt Hamer et al., 2005) - [PNAS Full Text](https://www.pnas.org/doi/10.1073/pnas.0502887102).

### 2. Therapy-Driven Evolution (LITT)
Therapies such as radiation and thermal ablation do not just reduce tumor bulk; they drive selection.
* **LITT Mechanism:** LITT utilizes a laser fiber to deliver thermal energy, creating a central zone of necrosis ($>60^\circ$C) surrounded by a sublethal thermal penumbra ($43-60^\circ$C).
* **Penumbra Effects:** Surviving cells in this penumbra experience heat stress, transient blood-brain barrier (BBB) disruption, and hypoxia, likely driving distinct transcriptional states and eventual recurrence.
* **Model Validation:** The Nagaraja et al. model confirms that while LITT achieves near-complete ablation of the central mass, viable tumor cells persist in the periphery, serving as the seed for recurrence.

### 3. Primary vs. Recurrent GBM Signatures
Matched patient datasets provide a template for analyzing recurrence.
* **Recurrence Signatures:** Recurrent GBM consistently shows upregulation of mesenchymal/stromal programs, myelination, and immune interactions (e.g., Fcγ receptor, complement).
* **Downregulation:** Purely proliferative and cell-cycle pathways are often downregulated in recurrence compared to primary tumors.

**Key Reference:**
* *Multidimensional analysis of matched primary and recurrent glioblastoma...* (2025) - [J Neuropathol Exp Neurol](https://academic.oup.com/jnen/article/84/1/45/7826743).

---

## The Gap: Why This Dataset is Needed

Despite the establishment of LITT as a therapy and U251 as a model, current literature lacks a unified transcriptomic study that:
1.  **Benchmarks** U251 culture-to-brain adaptation at RNAseq resolution.
2.  **Profiles** the orthotopic tumor specifically after LITT focal ablation.
3.  **Longitudinally samples** the recurrence to test for convergence on mesenchymal/stress-tolerant states.

This repository houses the analysis and data to address this unmet need, integrating differential gene expression and pathway analysis across all three evolutionary stages.

---

## Bioinformatics Workflow

To ensure reproducibility and standardization, this project utilizes the [nf-core/rnaseq](https://nf-co.re/rnaseq) analysis pipeline.

### Pipeline Steps
The workflow processes raw FASTQ files through the following steps:
1.  **QC:** FastQC and MultiQC.
2.  **Trimming:** TrimGalore! / Cutadapt.
3.  **Alignment:** STAR aligner mapped to the human reference genome (GRCh38).
4.  **Quantification:** Salmon / RSEM for gene-level counts.
5.  **Differential Expression:** Downstream analysis performed in R (DESeq2) comparing:
    * *In Vitro* (C2B) vs. *Primary* (IL64B, IL67B, IL68B, IL69B) (Adaptation signature)
    * *Primary* vs. *Recurrent* (IL66B, NL70B, NL71B) (Therapy-driven evolution)

### Usage
To reproduce the processing of raw data:

```bash
nextflow run nf-core/rnaseq \
    -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --aligner star_salmon
