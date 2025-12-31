#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH -o preview_%j.out
#SBATCH -e preview_%j.err

source "${HOME}/mambaforge/etc/profile.d/mamba.sh"
source activate nextflow

export NXF_SINGULARITY_CACHEDIR=${HOME}/singularity_cache
export XDG_RUNTIME_DIR=${HOME}/xdr
export NXF_EXECUTOR=slurm
export NXF_OPTS="-Xms2G -Xmx8G"

# Using -preview flag
nextflow run nf-core/rnaseq \
    -r 3.14.2 \
    -preview \
    -profile singularity,slurm \
    --input ANALYSIS/samplesheet.csv \
    --outdir ANALYSIS/results \
    --genome GRCh38 \
    --bbsplit_fasta_list ANALYSIS/bbsplit.csv \
    --skip_bbsplit false \
    --save_bbsplit_reads \
    --max_cpus 16 \
    --max_memory '64.GB'
