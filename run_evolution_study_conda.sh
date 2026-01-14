#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=go2432@wayne.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o xengsort_prod_%j.out
#SBATCH -e xengsort_prod_%j.err

# Trace every command to the log file for debugging
set -x

# ==========================================
# 1. ENVIRONMENT SETUP (Conda Activate Method)
# ==========================================
echo "LOG: Sourcing and Activating Conda..."
source "${HOME}/mambaforge/etc/profile.d/conda.sh"
conda activate nextflow

# Sanity Check
echo "LOG: Verifying Nextflow..."
which nextflow
nextflow -v

# ==========================================
# 2. CACHE & STORAGE SETUP
# ==========================================
# Use local home directory for small config/cache files
export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
mkdir -p $XDG_RUNTIME_DIR $NXF_SINGULARITY_CACHEDIR

# Use SCRATCH for the heavy lifting (intermediate files)
WORK_DIR="$(pwd)/work" 

# ==========================================
# 3. SAFETY LOCKS (Crucial for RSeQC)
# ==========================================
# These variables leak from Conda into Singularity and crash R/Python tools.
# Unsetting them forces the container to use its own internal libraries.
export NXF_SINGULARITY_HOME_MOUNT=true
unset LD_LIBRARY_PATH
unset PYTHONPATH
unset R_LIBS
unset R_LIBS_USER
unset R_LIBS_SITE

# ==========================================
# 4. PIPELINE EXECUTION
# ==========================================

# --- STEP 1: XENGSORT (Deconvolution) ---
echo "RUNNING STEP 1: XENGSORT"
nextflow run main.nf -profile singularity \
    --input "ANALYSIS/samplesheet.csv" \
    --host_fasta "/wsu/home/go/go24/go2432/u251-transcriptomic-evolution/ANALYSIS/refs/rat/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz" \
    --graft_fasta "/wsu/home/go/go24/go2432/u251-transcriptomic-evolution/ANALYSIS/refs/human/GRCh38.primary_assembly.genome.fa.gz" \
    --outdir "ANALYSIS/xengsort_out" \
    -w "${WORK_DIR}" \
    -resume -ansi-log false

unset NXF_PARAMS

# --- STEP 2: RNA-SEQ (Quantification & QC) ---
echo "RUNNING STEP 2: RNA-SEQ"
nextflow run nf-core/rnaseq \
    -r 3.14.0 \
    -profile singularity \
    --input "$(pwd)/ANALYSIS/human_evolution_paired.csv" \
    --outdir "$(pwd)/ANALYSIS/results_human_final" \
    --fasta "$(pwd)/ANALYSIS/refs/human/GRCh38.primary_assembly.genome.fa.gz" \
    --gtf "$(pwd)/ANALYSIS/refs/human/GRCh38.primary_assembly.annotation.gtf.gz" \
    --skip_bbsplit true \
    --save_unaligned \
    --aligner star_salmon \
    --max_cpus 16 \
    --max_memory '120.GB' \
    -w "${WORK_DIR}" \
    -resume \
    -ansi-log false

unset NXF_PARAMS

# --- STEP 3: DIFFERENTIAL ABUNDANCE ---
echo "RUNNING STEP 3: DIFFERENTIAL ABUNDANCE"
nextflow run nf-core/differentialabundance \
    -r 1.5.0 \
    -profile singularity \
    --input "$(pwd)/ANALYSIS/metadata_therapy.csv" \
    --contrasts "$(pwd)/ANALYSIS/contrasts_therapy.csv" \
    --matrix "$(pwd)/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv" \
    --transcript_length_matrix "$(pwd)/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_lengths.tsv" \
    --gtf "$(pwd)/ANALYSIS/refs/human/GRCh38.primary_assembly.annotation.gtf.gz" \
    --exploratory_main_variable Classification \
    --shinyngs_build_app \
    -params-file therapy_params.yaml \
    -w "${WORK_DIR}" \
    -resume \
    -ansi-log false
