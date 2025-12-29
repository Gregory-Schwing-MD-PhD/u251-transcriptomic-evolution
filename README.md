# Transcriptomic Evolution of U251 Glioblastoma Cells: From Culture to Orthotopic Growth and Post-therapy Recurrence

![Visual Abstract](./ASSETS/Visual_abstract.jpg)

## Project Overview

This project addresses a critical gap in glioblastoma (GBM) modeling: the lack of longitudinal RNA-sequencing (RNAseq) datasets tracking tumor evolution through clinically relevant stages. [cite_start]Specifically, we investigate the transcriptomic shifts of U251 cells from **standard in vitro culture**, to **orthotopic growth** in the murine brain, and finally through **recurrence following Laser Interstitial Thermal Therapy (LITT)**[cite: 1, 4, 49].

[cite_start]While U251 is a widely used model, classic profiling shows that gene expression shifts substantially when implanted orthotopically, underscoring the dominance of the brain microenvironment[cite: 2, 3]. This project aims to quantify these adaptations and define the transcriptional consequences of focal thermal ablation.

---

## Background & Scientific Context

### 1. In Vitro vs. Orthotopic Xenografts
[cite_start]Evidence suggests that the brain microenvironment exerts a dominant influence on transcriptional states[cite: 8].
* [cite_start]**Convergent Profiles:** Gene expression profiles of GBM lines (U251, U87) become more similar to each other—and to patient GBM—in orthotopic settings than in vitro[cite: 8, 11].
* [cite_start]**Key Pathways:** Orthotopic tumors upregulate extracellular matrix (ECM) remodeling, cell adhesion, and angiogenesis, whereas in vitro cultures favor classical proliferation genes[cite: 9].
* [cite_start]**Clinical Relevance:** Orthotopic xenografts cluster closely with patient GBM samples, recapitulating hypoxia and invasion signatures that subcutaneous models fail to capture[cite: 11].

**Key References:**
* [cite_start]*Influence of in vivo growth on human glioma cell line gene expression* (De Witt Hamer et al., 2005) - [PNAS Full Text](https://www.pnas.org/doi/10.1073/pnas.0502887102)[cite: 6].
* [cite_start]*Molecular profiling indicates orthotopic xenograft... is a more clinically relevant model* (2012) - [Clin Cancer Res Full Text](https://pmc.ncbi.nlm.nih.gov/articles/PMC3164941/)[cite: 10].

### 2. Therapy-Driven Evolution (Radiation & LITT)
[cite_start]Therapies such as radiation do not just reduce tumor bulk; they drive selection[cite: 17].
* [cite_start]**Radiation Effects:** Fractionated radiation selects for subclones with mesenchymal, inflammatory, and stress-tolerant features[cite: 15, 16].
* **LITT Mechanism:** LITT creates a central zone of necrosis and a sublethal thermal penumbra. [cite_start]Surviving cells in this penumbra experience heat stress and hypoxia, likely driving distinct transcriptional states[cite: 24, 25].
* [cite_start]**Hypothesis:** Recurrent tumors after LITT are expected to enrich for mesenchymal and DNA-damage-response signatures, analogous to radiation-induced evolution[cite: 18].

**Key References:**
* [cite_start]*Radiation Drives the Evolution of Orthotopic Xenografts...* (2019) - [Cancer Res Full Text](https://aacrjournals.org/cancerres/article/79/23/6032/640077/Radiation-Drives-the-Evolution-of-Orthotopic)[cite: 15].
* [cite_start]*Laser Interstitial Thermal Therapy for Recurrent Glioblastoma* (2021) - [Neurosurgery Full Text](https://pmc.ncbi.nlm.nih.gov/articles/PMC11605670/)[cite: 22].

### 3. Primary vs. Recurrent GBM Signatures
[cite_start]Matched patient datasets provide a template for analyzing recurrence[cite: 29].
* [cite_start]**Recurrence Signatures:** Recurrent GBM consistently shows upregulation of mesenchymal/stromal programs, myelination, and immune interactions (e.g., Fcγ receptor, complement)[cite: 31, 32, 33].
* [cite_start]**Downregulation:** Purely proliferative and cell-cycle pathways are often downregulated in recurrence compared to primary tumors[cite: 31, 35].

**Key Reference:**
* [cite_start]*Multidimensional analysis of matched primary and recurrent glioblastoma...* (2025) - [J Neuropathol Exp Neurol](https://academic.oup.com/jnen/article/84/1/45/7826743)[cite: 30].

---

## The Gap: Why This Dataset is Needed

Despite the establishment of LITT as a therapy and U251 as a model, current literature lacks a unified transcriptomic study that:
1.  [cite_start]**Benchmarks** U251 culture-to-brain adaptation at RNAseq resolution[cite: 4, 49].
2.  [cite_start]**Profiles** the orthotopic tumor specifically after LITT focal ablation[cite: 26].
3.  [cite_start]**Longitudinally samples** the recurrence to test for convergence on mesenchymal/stress-tolerant states[cite: 50].

[cite_start]This repository houses the analysis and data to address this unmet need, integrating differential gene expression and pathway analysis across all three evolutionary stages[cite: 51].

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
    * *In Vitro* vs. *Orthotopic Pre-LITT* (Adaptation signature)
    * *Orthotopic Pre-LITT* vs. *Post-LITT Recurrence* (Therapy-driven evolution)

### Usage
To reproduce the processing of raw data:

```bash
nextflow run nf-core/rnaseq \
    -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --aligner star_salmon
