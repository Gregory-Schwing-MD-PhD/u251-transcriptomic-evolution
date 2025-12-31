#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2      # Only 2 needed; it just manages other jobs
#SBATCH --mem=12G              # Matches your NXF_OPTS overhead
#SBATCH --time=7-00:00:00      # Keep it long; this "Master" job must stay alive
#SBATCH --mail-user=go2432@wayne.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o rnaseq_master_%j.out
#SBATCH -e rnaseq_master_%j.err

# 1. Environment Setup
source "${HOME}/mambaforge/etc/profile.d/mamba.sh"
source activate nextflow

# 2. Singularity & Path Setup
export NXF_SINGULARITY_CACHEDIR=${HOME}/singularity_cache
mkdir -p $NXF_SINGULARITY_CACHEDIR
mkdir -p ${HOME}/xdr
export XDG_RUNTIME_DIR=${HOME}/xdr

# 3. Nextflow Tuning (from your successful runs)
export NXF_EXECUTOR=slurm
export NXF_OPTS="-Xms2G -Xmx8G"

# 4. The Command
# Note: We use -profile singularity,slurm to combine tool handling + job submission
nextflow run nf-core/rnaseq \
    -r 3.14.2 \
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
