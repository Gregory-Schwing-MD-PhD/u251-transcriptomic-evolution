#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --mail-user=go2432@wayne.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o de_prod_%j.out
#SBATCH -e de_prod_%j.err

# 1. Environment Setup
source "${HOME}/mambaforge/etc/profile.d/conda.sh"
source activate nextflow
export PATH="${HOME}/mambaforge/envs/nextflow/bin:$PATH"
unset JAVA_HOME

# 2. Path Setup
export XDG_RUNTIME_DIR="${HOME}/xdr"
mkdir -p $XDG_RUNTIME_DIR

# 3. The Production Command
nextflow run nf-core/differentialabundance \
    -r 1.5.0 \
    -profile singularity \
    --input "$(pwd)/ANALYSIS/metadata.csv" \
    --contrasts "$(pwd)/ANALYSIS/contrasts.csv" \
    --matrix "$(pwd)/ANALYSIS/results/star_salmon/salmon.merged.gene_counts.tsv" \
    --transcript_length_matrix "$(pwd)/ANALYSIS/results/star_salmon/salmon.merged.gene_lengths.tsv" \
    --genome GRCh38 \
    --gene_sets_files "$(pwd)/ANALYSIS/refs/h.all.v2023.2.Hs.symbols.gmt" \
    --study_name "U251_LITT_Evolution" \
    --outdir "$(pwd)/ANALYSIS/results_differential" \
    --shinyngs_build_app \
    --deseq2_min_replicates_for_replace 3 \
    --gsea_run true \
    --gprofiler2_run true \
    --gprofiler2_organism hsapiens \
    -resume
