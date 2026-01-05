#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=go2432@wayne.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o rnaseq_prod_%j.out
#SBATCH -e rnaseq_prod_%j.err

# 1. Environment Setup
source "${HOME}/mambaforge/etc/profile.d/conda.sh"
source activate nextflow
export PATH="${HOME}/mambaforge/envs/nextflow/bin:$PATH"
unset JAVA_HOME

# 2. Path Setup
export XDG_RUNTIME_DIR="${HOME}/xdr"
mkdir -p $XDG_RUNTIME_DIR

# 3. The Production Command
# Note: '-profile slurm' is removed to avoid the "Unknown profile" error.
# Note: config file handles executor, singularity, and resource limits.
nextflow run nf-core/rnaseq \
    -r 3.22.2 \
    -profile singularity \
    --input "$(pwd)/ANALYSIS/samplesheet.csv" \
    --outdir "$(pwd)/ANALYSIS/results" \
    --fasta "$(pwd)/ANALYSIS/refs/human/GRCh38.primary_assembly.genome.fa.gz" \
    --gtf "$(pwd)/ANALYSIS/refs/human/GRCh38.primary_assembly.annotation.gtf.gz" \
    --bbsplit_fasta_list "$(pwd)/ANALYSIS/bbsplit.csv" \
    --skip_bbsplit false \
    --save_bbsplit_reads \
    --max_cpus 16 \
    --max_memory '64.GB' \
    -resume
