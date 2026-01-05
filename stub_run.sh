#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH -o stub_%j.out
#SBATCH -e stub_%j.err

# 1. Environment Setup (The Verified Way)
source "${HOME}/mambaforge/etc/profile.d/conda.sh"
source activate nextflow
export PATH="${HOME}/mambaforge/envs/nextflow/bin:$PATH"
unset JAVA_HOME

# 2. Paths
# XDG is important for Singularity stability on clusters
export XDG_RUNTIME_DIR="${HOME}/xdr"
mkdir -p $XDG_RUNTIME_DIR

# 3. Stub-run Command
# Nextflow picks up the executor and QoS from 'nextflow.config' automatically
nextflow run nf-core/rnaseq \
    -r 3.22.2 \
    -stub-run \
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
