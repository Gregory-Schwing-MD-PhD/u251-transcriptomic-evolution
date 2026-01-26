# Reference Dataset: Neftel et al. (2019)

This directory contains the single-cell RNA-seq reference data used to generate the "Signature Matrix" for CIBERSORTx deconvolution.

## Source Information
* **Study:** An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma
* **Authors:** Neftel C, Laffy J, Filbin MG, et al.
* **Journal:** *Cell*, 2019
* **DOI:** [10.1016/j.cell.2019.06.055](https://doi.org/10.1016/j.cell.2019.06.055)
* **Accession:** [Broad Single Cell Portal (SCP393)](https://singlecell.broadinstitute.org/single_cell/study/SCP393/single-cell-rna-seq-of-adult-and-pediatric-glioblastoma)

## Files Description

| File Name | Description | Role in Pipeline |
| :--- | :--- | :--- |
| **`IDHwtGBM.processed.SS2.logTPM.txt.gz`** | Processed expression matrix (log2 TPM) from Smart-Seq2 full-length protocol. Filtered for IDH-wildtype Glioblastoma. | **Input Expression:** Converted to linear TPM by `prepare_cibersort.R` to build the signature matrix. |
| **`IDHwt.GBM.Metadata.SS2.txt`** | Clinical and cellular metadata, including the 4-state classification scores (MES, AC, OPC, NPC). | **Input Phenotypes:** Used to assign "Truth" labels to cells for training the CIBERSORTx signature. |
| **`.gitkeep`** | Empty invisible file. | **Git Tracking:** Ensures this directory structure persists in the repo even if data files are ignored. |

## Download Instructions (Abbreviated)

**Note:** Direct download links are not available due to access control. You must log in.

1.  **Navigate:** Go to the [Broad Single Cell Portal Study SCP393](https://singlecell.broadinstitute.org/single_cell/study/SCP393/single-cell-rna-seq-of-adult-and-pediatric-glioblastoma).
2.  **Login:** Sign in using a Google account (required for download).
3.  **Select Tab:** Click the **"Download"** tab.
4.  **Download Files:**
    * Find the **Expression** section and download `IDHwtGBM.processed.SS2.logTPM.txt.gz`.
    * Find the **Metadata** section and download `IDHwt.GBM.Metadata.SS2.txt`.
5.  **Placement:** Move both files into `ANALYSIS/refs/neftel_2019/`.

## Usage Notes
* **Do not unzip manually.** The pipeline script (`prepare_cibersort.R`) is designed to read the `.gz` file directly.
* **Data Privacy:** These files are derived from human patients. While de-identified, ensure compliance with data handling protocols.

## Pipeline Integration
These files are automatically detected and processed by `cibersortx.nf` during **Step 3.5** of the evolution study.
