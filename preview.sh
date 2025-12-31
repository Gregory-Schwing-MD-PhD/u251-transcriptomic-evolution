#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH -o preview_%j.out
#SBATCH -e preview_%j.err

# 1. Environment Setup (The Verified Way)
source "${HOME}/mambaforge/etc/profile.d/conda.sh"
source activate nextflow
export PATH="${HOME}/mambaforge/envs/nextflow/bin:$PATH"
unset JAVA_HOME

# 2. Paths
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_EXECUTOR=slurm

# 3. Preview Command
nextflow run nf-core/rnaseq \
    -r 3.22.2 \
    -profile singularity \
    --input ANALYSIS/samplesheet.csv \
    --outdir ANALYSIS/results \
    --genome GRCh38 \
    --bbsplit_fasta_list ANALYSIS/bbsplit.csv \
    --skip_bbsplit false \
    --save_bbsplit_reads \
    --max_cpus 16 \
    --max_memory '64.GB' \
    -resume
