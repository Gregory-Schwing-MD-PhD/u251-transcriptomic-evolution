#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2      # Master job only needs 2 CPUs to orchestrate
#SBATCH --mem=16G              # Master job overhead
#SBATCH --time=7-00:00:00      # 7-day limit (Master must stay alive while workers run)
#SBATCH --mail-user=go2432@wayne.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o rnaseq_prod_%j.out
#SBATCH -e rnaseq_prod_%j.err

# 1. Environment Setup (The Verified "Quiet" Way)
source "${HOME}/mambaforge/etc/profile.d/conda.sh"
source activate nextflow
export PATH="${HOME}/mambaforge/envs/nextflow/bin:$PATH"
unset JAVA_HOME

# 2. Paths & Cluster Settings
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_EXECUTOR=slurm

# Ensure directories exist
mkdir -p $NXF_SINGULARITY_CACHEDIR
mkdir -p $XDG_RUNTIME_DIR

# 3. The Production Command
# Using -resume is critical for large runs
nextflow run nf-core/rnaseq \
    -r 3.22.2 \
    -profile singularity,slurm \
    --input ANALYSIS/samplesheet.csv \
    --outdir ANALYSIS/results \
    --genome GRCh38 \
    --bbsplit_fasta_list ANALYSIS/bbsplit.csv \
    --skip_bbsplit false \
    --save_bbsplit_reads \
    --max_cpus 16 \
    --max_memory '64.GB' \
    -resume
