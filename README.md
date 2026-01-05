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
* [cite_start]**LITT Mechanism:** LITT utilizes a laser fiber to deliver thermal energy, creating a central zone of necrosis ($>60^\circ$C) surrounded by a sublethal thermal penumbra ($43-60^\circ$C)[cite: 63, 65, 68].
* [cite_start]**Penumbra Effects:** Surviving cells in this penumbra experience heat stress, transient blood-brain barrier (BBB) disruption, and hypoxia, likely driving distinct transcriptional states and eventual recurrence[cite: 65, 210, 214].
* [cite_start]**Model Validation:** The Nagaraja et al. model confirms that while LITT achieves near-complete ablation of the central mass, viable tumor cells persist in the periphery, serving as the seed for recurrence[cite: 210, 235].

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
4.  **Simultaneously profiles** the host microenvironment response (stroma/immune) to distinguish tumor-intrinsic evolution from extrinsic tissue scarring.

This repository houses the analysis and data to address this unmet need, integrating differential gene expression and pathway analysis across all three evolutionary stages.

---

## Bioinformatics Workflow: The Dual-Species Strategy

To maximize the utility of the orthotopic xenograft model, this project employs a **Dual-Species Workflow**. Raw sequencing reads are computationally sorted into **Human (Tumor)** and **Rat (Host)** streams, creating two parallel experiments from a single dataset.

| Experiment | **1. Tumor Evolution** | **2. Host Microenvironment** |
| :--- | :--- | :--- |
| **Target Organism** | Human (U251 Cells) | Rat (Brain Stroma/Microglia) |
| **Input Data** | Human-mapping reads | Rat-mapping reads (discarded from Exp 1) |
| **Biological Goal** | Track tumor adaptation & resistance | Track inflammation, gliosis & scarring |
| **Key Contrast** | Primary vs. Recurrent (Resistance) | Tumor vs. Control (Inflammation) |

### Phase 1: Alignment & Species Sorting
**Tool:** `nf-core/rnaseq` (v3.22.2)

A critical challenge in orthotopic xenografts is the high sequence conservation between human tumor cells and the rat host brain. To address this, we employ **Competitive Alignment** using [BBSplit](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbsplit-guide/).
1.  **Dual Indexing:** Reads are mapped simultaneously to a combined reference of **Human (GRCh38)** and **Rat (mRatBN7.2)**.
2.  **Disambiguation:** Reads mapping best to Rat are segregated; reads mapping best to Human are retained.
3.  **Quantification:** "Clean" human reads are aligned (STAR) and quantified (Salmon).

*Note: Reference genomes are loaded locally from `ANALYSIS/refs/` to ensure reproducibility and prevent version drift.*

```bash
nextflow run nf-core/rnaseq \
    -r 3.22.2 \
    -profile singularity \
    --input ANALYSIS/samplesheet.csv \
    --outdir ANALYSIS/results \
    --fasta ANALYSIS/refs/human/GRCh38.primary_assembly.genome.fa \
    --gtf ANALYSIS/refs/human/GRCh38.primary_assembly.annotation.gtf \
    --bbsplit_fasta_list ANALYSIS/bbsplit.csv \
    --skip_bbsplit false \
    --save_bbsplit_reads \
    --max_cpus 16 \
    --max_memory '64.GB' \
    -resume
```

### Phase 2A: Tumor Analysis (Human Reads)
**Tool:** `nf-core/differentialabundance` (v1.5.0)

Raw counts are transformed into biological insights using DESeq2. This phase applies generalized linear modeling to identify genes that are significantly perturbed by the microenvironment or therapeutic stress.

**Rationale for Contrasts**
1.  **Brain Adaptation** (*Culture vs. Primary*): Defines the "Engraftment Shock." Identifies genes required to transition from plastic to the brain parenchyma.
2.  **Therapy Impact** (*Primary vs. Recurrent*): Defines "Resistance." Isolates the specific transcriptomic shifts driven by thermal ablation and recovery.
3.  **Core Brain Signature** (*In Vitro vs. All In Vivo*): Defines "Invasion." By grouping Primary and Recurrent tumors against Culture, we identify the universal machinery required for U251 survival in the brain, independent of therapy.

```bash
nextflow run nf-core/differentialabundance \
    -r 1.5.0 \
    -profile singularity \
    --input "$(pwd)/ANALYSIS/metadata.csv" \
    --contrasts "$(pwd)/ANALYSIS/contrasts.csv" \
    --matrix "$(pwd)/ANALYSIS/results/star_salmon/salmon.merged.gene_counts.tsv" \
    --transcript_length_matrix "$(pwd)/ANALYSIS/results/star_salmon/salmon.merged.gene_lengths.tsv" \
    --gtf "$(pwd)/ANALYSIS/refs/human/GRCh38.primary_assembly.annotation.gtf" \
    --genesets "$(pwd)/ANALYSIS/refs/h.all.v2023.2.Hs.symbols.gmt" \
    --study_name "U251_LITT_Evolution" \
    --outdir "$(pwd)/ANALYSIS/results_differential" \
    --shinyngs_build_app \
    --deseq2_min_replicates_for_replace 3 \
    -c gsea_fix.config \
    -resume
```

### Phase 2B: Host Microenvironment Analysis (Rat Reads)
**Tool:** `nf-core/rnaseq` (Re-run) & `nf-core/differentialabundance`

Reads identified as "Rat" by BBSplit in Phase 1 are not discarded but are re-analyzed to profile the host response. This "inverted" analysis treats the tumor samples as "Inflamed Brain" and the control samples as "Healthy Brain."

**Rationale for Contrasts**
1.  **Tumor Inflammation** (*Primary Tumor vs. Control Brain*): Identifies microglial activation (*Aif1*), astrogliosis (*Gfap*), and vascular recruitment driven by the tumor presence.
2.  **Ablation Scarring** (*Recurrent Tumor vs. Primary Tumor*): Differentiates the chronic inflammatory signature of the LITT burn scar from the baseline tumor inflammation.

**Step 1: Re-align Rat Reads**
*Note: We map the discarded reads to the local Rat reference (`mRatBN7.2`) using the same resource constraints as the primary analysis.*

```bash
nextflow run nf-core/rnaseq \
    -r 3.22.2 \
    -profile singularity \
    --input ANALYSIS/samplesheet_host.csv \
    --outdir ANALYSIS/results_host \
    --fasta ANALYSIS/refs/rat/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa \
    --gtf ANALYSIS/refs/rat/Rattus_norvegicus.mRatBN7.2.110.gtf \
    --remove_ribo_rna \
    --skip_bbsplit \
    --max_cpus 16 \
    --max_memory '64.GB' \
    -resume
```

**Step 2: Differential Abundance (Host)**
*Note: Requires the Rat-specific Hallmark GMT file generated via `msigdbr`.*

```bash
nextflow run nf-core/differentialabundance \
    -r 1.5.0 \
    -profile singularity \
    --input "$(pwd)/ANALYSIS/metadata_host.csv" \
    --contrasts "$(pwd)/ANALYSIS/contrasts_host.csv" \
    --matrix "$(pwd)/ANALYSIS/results_host/star_salmon/salmon.merged.gene_counts.tsv" \
    --transcript_length_matrix "$(pwd)/ANALYSIS/results_host/star_salmon/salmon.merged.gene_lengths.tsv" \
    --gtf "$(pwd)/ANALYSIS/refs/rat/Rattus_norvegicus.mRatBN7.2.110.gtf" \
    --genesets "$(pwd)/ANALYSIS/refs/r.hallmark.v2023.2.Rn.symbols.gmt" \
    --study_name "U251_Host_Response" \
    --outdir "$(pwd)/ANALYSIS/results_host_differential" \
    --shinyngs_build_app \
    --deseq2_min_replicates_for_replace 3
```

## 🛠 Software Requirements & Installation

This pipeline uses a containerized infrastructure to ensure reproducibility. The core workflow is managed by **Nextflow**, software dependencies are isolated via **Singularity**, and the local runtime environment is managed by **Mamba**.

### Core Technologies
* **[Nextflow](https://www.nextflow.io/) (v23.10.0):** Orchestrates data flow, manages Slurm job submissions, and handles error recovery via the `-resume` flag.
* **[Singularity/Apptainer](https://apptainer.org/) (>=v3.6):** Executes the bioinformatics tools (STAR, Salmon, BBSplit) within isolated containers to ensure version consistency across cluster nodes.
* **[Mamba](https://mamba.readthedocs.io/):** A fast implementation of Conda used to manage the Nextflow installation and environment.
* **[Bioconda](https://bioconda.github.io/):** The software channel providing the bioinformatics-specific packages.

---

### 🚀 Setup and Installation

Run the following commands to install the package manager, build the environment, and activate it:

# 1. Install Mambaforge (if not already present)
```bash
# Download the installer
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"

# Run the installer
bash Mambaforge-$(uname)-$(uname -m).sh

# Refresh your profile
source ~/.bashrc
```

# 2. Create the environment from the yml file
```bash
mamba env create -f envs/nextflow.yml
```

# 3. Activate the environment
```bash
mamba activate nextflow
```

# 4. Verify installation
```bash
nextflow -v
singularity --version
```

## Additional Resources

### Reference Genomes & Databases
For alignment and enrichment analysis, the following specific assemblies and databases were used. These must be downloaded locally to `ANALYSIS/refs/`.

**1. Human Reference (Target - Tumor)**
* **Assembly:** GRCh38 (GENCODE Release 44)
* **Sequence (FASTA):**
    * [Download (ftp.ebi.ac.uk)](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz)
* **Annotation (GTF):**
    * [Download (ftp.ebi.ac.uk)](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.annotation.gtf.gz)

**2. Rat Reference (Host - Microenvironment)**
* **Assembly:** mRatBN7.2 (Ensembl Release 110)
* **Sequence (FASTA):**
    * [Download (ftp.ensembl.org)](http://ftp.ensembl.org/pub/release-110/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz)
* **Annotation (GTF):**
    * [Download (ftp.ensembl.org)](http://ftp.ensembl.org/pub/release-110/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.110.gtf.gz)

**3. Pathway Databases**
* **Human Hallmarks (Standard):**
    * [Download MSigDB v2023.2 (Broad Institute)](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt)
* **Rat Hallmarks (Generated):**
    * *Note:* There is no official Rat GMT. The file `r.hallmark.v2023.2.Rn.symbols.gmt` must be generated locally using the R package `msigdbr` to map Human Hallmarks to Rat orthologs.
* **Brain.GMT (Tissue-Specific):**
    * A curated database of brain-related functional gene sets (e.g., Glioblastoma signaling, Microglial activation, Astrocyte polarization) not found in standard MSigDB.
    * *Source:* Pusarla et al., 2024 ("A Curated Database of Brain-Related Functional Gene Sets").
    * [Access Database / Paper](https://pubmed.ncbi.nlm.nih.gov/39049932/)

### Tools & Documentation
* [nf-core/rnaseq Documentation](https://nf-co.re/rnaseq)
* [nf-core/differentialabundance Documentation](https://nf-co.re/differentialabundance/1.5.0)
* [BBTools / BBSplit Guide](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbsplit-guide/)
* [Visualase® Clinical System Info](https://www.medtronic.com/us-en/healthcare-professionals/products/neurological/laser-ablation-systems/visualase.html)
